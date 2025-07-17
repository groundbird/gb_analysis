#!/usr/bin/env python
import datetime
import django
import sys
import os
import json
import scipy
from scipy import signal
import numpy as np

from pathlib import Path
from datetime import timezone, timedelta
from stat import filemode
from os import getuid, getgroups, getgid, getpid
from argparse import ArgumentParser
from copy import deepcopy

from kidslist import KidsList
from mkid_pylibs.readfile import readfile_swp

DATA_DIR = Path('/Users/pika/GroundBIRD/data')

# Django preparation
MEASWEB_PATH = Path('/Users/pika/GroundBIRD/meas_web')
sys.path.append(str(MEASWEB_PATH))
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "meas_web.settings")
django.setup()

from file_manager.models import Tod, Sweep, RheaTod, RheaSweep, Measurement

# Signal generator
SG_PATH = Path('/Users/pika/GroundBIRD/quicksyn')
sys.path.append(str(SG_PATH))
from sg_manager import QuickSyn

# rhea_comm
RHEA_PATH = Path('/Users/pika/GroundBIRD/analyzer_db/rhea_comm')
sys.path.append(str(RHEA_PATH))
#from fpga_control import fpga_control # for rhea_comm_old_202209
from fpga_control import FPGAControl
from measure_mulswp import measure_mulswp
from measure_tod import measure_tod
#from tod_info import TodInfo, TodInfoException

def two_div(val):
    return len(bin(val-1)) - 2

def is_writable(path):
    p_st = path.stat()
    fm_str = filemode(p_st.st_mode)

    if fm_str[2] == 'w':# Write access to the owner
        uid = getuid()
        if uid == p_st.st_uid:
            return True

    if fm_str[5] == 'w':# Write access to the group
        if p_st.st_gid in getgroups():
            return True

    if fm_str[8] == 'w':# Write access to everyone
        return True

    return False

def path_checker(path):
    if not path.exists():
        raise RuntimeError(f'Path {path} does not exist.')

    if not path.is_dir():
        raise RuntimeError(f'Path {path} is not a directory.')

    if not is_writable(path):
        raise RuntimeError(f'You do not have a write access to the path {path}')

def path_creator(dirpath):
    utcnow = datetime.datetime.now(tz=timezone.utc)
    d = dirpath.joinpath(f'{utcnow.year:04d}{utcnow.month:02d}{utcnow.day:02d}')
    d = d.joinpath(f'data_{utcnow.hour:02d}{utcnow.minute:02d}{utcnow.second:02d}')
    d.mkdir(exist_ok=True, parents=True)
    return d

def get_int_or_float(v):
    number_as_float = float(v)
    number_as_int = int(number_as_float)
    return number_as_int if number_as_float == number_as_int else number_as_float

def make_config(fn):
    import configparser
    ini = configparser.ConfigParser()
    ini.optionxform = str
    ini.read(fn)

    ret = {}
    for isect in ini.sections():
        confitems = config_cleaner(ini.items(isect))
        if 'config' in confitems.keys():
            try:
                subini = configparser.ConfigParser()
                subini.optionxform = str
                subini.read(confitems['config'])
                subconfitems = config_cleaner(subini.items(isect))
                subconfitems.update(confitems)
                ret[isect] = subconfitems
            except:
                print("Failed to read \"{isect}\" in sub configuration file: {confitems['config']}")
                ret[isect] =confitems
                pass
            pass
        else:
            ret[isect] =confitems
        pass
    return ret

def config_cleaner(items):
    newitems = {}
    for (k,v) in items:
        v = v.strip()
        if (v[0] == '"' or v[0] == "'") and v[-1] == v[0]:
            v = v[1:-1]

        if v.upper() == 'TRUE':
            ret = True
        elif v.upper() == 'FALSE':
            ret = False
        elif v.upper() == 'NONE':
            ret = None
        else:
            try:
                if '.' in v:
                    ret = float(v)
                else:
                    ret = int(v)
            except:
                ret = v
        newitems[k] = ret
    return newitems

def main(klist_path, length,
         fpga=None,
         swap_dac=False,
         swap_adc=False,
         width_sweep=3.0,
         step_sweep=0.01,
         path_dir=None,
         sg_config={},
         rate_kHz=1,
         comment='Created by meas_db.py',
         nodb=False,
         label=''):

    klist = KidsList(klist_path)

#    MAX_CH = fpga.MAX_CH # for rhea_comm_old_202209
    max_ch = fpga.max_ch    
    input_len = 2**two_div(len(klist.tone_freqs))

    # power check
    power = klist.kids_power[0]
    for p in klist.kids_power + klist.blinds_power:
        if p != power:
            raise RuntimeError('This script currently does not support power adjustment')

#    if power < 1 or power*input_len > MAX_CH:
#        raise RuntimeError(f'exceeding max # of channels = {input_len}*{power} > {MAX_CH}')
# for rhea_comm_old_202209
    if power < 1 or power*input_len > max_ch:
        raise RuntimeError(f'exceeding max # of channels = {input_len}*{power} > {max_ch}')

    ############################################################################# Directory creation
    path_dir = path_dir

    if path_dir is None:
        path_dir = path_creator(DATA_DIR)
    else:
        path_dir = Path(path_dir)
        path_checker(path_dir)

    swp_fname = path_dir.joinpath(f'swp_{label}.rawdata')
    tod_fname = path_dir.joinpath(f'tod_{label}.rawdata')
    klist_fname = path_dir.joinpath(f'kids_{label}.list')

    ############################################################################# SG setting
    print(f'SG setting: {sg_config.items()}')
    print(f'            freq = {klist.sg_freq/1e9} GHz')
    quicksyn = QuickSyn(**sg_config)
    quicksyn.set_rfout(on=False)
    quicksyn.set_freq_mHz(int(klist.sg_freq*1000))
    quicksyn.set_rfout(on=True)
    quicksyn.close()

    ############################################################################# Sweep
    ## Registration to the database
    if not nodb:
        meas = Measurement(array=klist.array_name, klpath=klist_fname, comment=comment)
        meas.save()

        swp_dt = datetime.datetime.now(tz=timezone.utc)
        swp = RheaSweep(path=swp_fname, mtime=swp_dt, ntone=len(klist.tone_freqs), 
                        length=int(width_sweep/step_sweep), 
                        resolution=step_sweep*1e6, 
                        span=width_sweep*1e6,
                        sg=klist.sg_freq)
        swp.save()
        meas.sweeps.add(swp)
        meas.save()

    ## Measurement
    print(klist.tone_freqs)
    measure_mulswp(fpga      = fpga,
#                   MAX_CH    = MAX_CH, # for rhea_comm_old_202209
                   max_ch    = max_ch,                   
#                   f_centers = klist.tone_freqs,
#                   dds_f_MHz = klist.tone_freqs, # for rhea_comm_old_202209
                   dds_f_megahz = klist.tone_freqs,
                   width     = width_sweep,
                   step      = step_sweep,
                   fname     = str(swp_fname),
                   power     = power,
                   swap_dac  = swap_dac,
                   swap_adc  = swap_adc,
                   amps      = klist.tone_amps,
                   phases    = klist.tone_phases)

    ############################################################################# Fitting
    kldict_new = deepcopy(klist._kldict)

    swp_data = readfile_swp('rhea', 
                            filename=str(swp_fname),
                            index=-1, # all
                            lo=klist.sg_freq)

    # 'blinds_freqs` will be added if only 'blinds_relfreqs' is assigned in kidslist
    is_newblinds = ('blinds_relfreqs' in kldict_new) and (not 'blinds_freqs' in kldict_new)
    if is_newblinds:
        kldict_new['blinds_freqs'] = [0]*len(klist.kids_index)

    def getfr_fit(swpd,fitconf):
        if not 'nfwhm' in fitconf: fitconf['nfwhm'] = -1
        if not 'frqrange' in fitconf:
            swpd.fitIQ(nfwhm=fitconf['nfwhm'])
        else:
            swpd.fitIQ(nfwhm=-1,frqrange=fitconf['frqrange'])
        return swpd.fitresult.params['fr'].value # RF, Hz

    def getfr_min(swpd):
        swp = swpd.amplitude
        f = swpd.f
        fcent = np.mean(f) # center frq
        swpmed = scipy.signal.savgol_filter(scipy.signal.medfilt(swp,3),5,1)
        is_minima = np.r_[True, swpmed[1:] < swpmed[:-1]] & np.r_[swpmed[:-1] < swpmed[1:], True] # find minima
        is_minima[0] = False  # not peak
        is_minima[-1] = False # not peak
        is_minima[is_minima] &= swp[is_minima] < swp.max() - swp.std()  # require `peak < (max - 1sigma)`
        ind_orig = np.array(np.arange(len(is_minima)).tolist()) # original index of freq
        ind = np.argmin(np.abs(fcent-f[is_minima])) # choose closest freq to center as fr
        ind_fr = ind_orig[is_minima][ind] # index of fr in original freq array
        if ind_fr>1 and ind_fr<len(swp)-2:
            ind_fr = np.argmin(swp[ind_fr-2:ind_fr+2]) + ind_fr-2 # remove filter effect
        return f[ind_fr] # RF, Hz
        
    for i, ki in enumerate(klist.kids_index):
        if klist.fit_conf[i] is False:
            fr = getfr_min(swp_data[ki])
        else:
            if klist.fit_conf[i] is True:
                fitconf = {'nfwhm':3}
            else:
                fitconf = klist.fit_conf[i]
            fr = getfr_fit(swp_data[ki],fitconf)

        fr_if = fr - klist.sg_freq # IF, Hz
        fr_if_kHz = int(fr_if/1e3 + 0.5) # IF, kHz
        kldict_new['kids_freqs'][i] = fr_if_kHz/1e3 # IF, MHz
        if is_newblinds:
            kldict_new['blinds_freqs'][i] = fr_if_kHz/1e3 + kldict_new['blinds_relfreqs'][i] # IF, MHz

    # for i, ki in enumerate(klist.kids_index):
    #     if klist.fit_conf[i] is False:
    #         print(f"KID#{i}: fit skipped.")
    #         continue
    #     elif klist.fit_conf[i] is True:
    #         fitconf = {'nfwhm':3}
    #     else:
    #         fitconf = klist.fit_conf[i]
    #         if not 'nfwhm' in fitconf: fitconf['nfwhm'] = -1
    #     swp = swp_data[ki].amplitude
    #     f = swp_data[ki].f
    #     fr = np.mean(f)
    #     a = scipy.signal.savgol_filter(scipy.signal.medfilt(swp,3),5,1)
    #     tmp = np.r_[True, a[1:] < a[:-1]] & np.r_[a[:-1] < a[1:], True]
    #     tmp[0] = False
    #     tmp[-1] = False
    #     tmp[tmp] &= swp[tmp] < swp.max() - swp.std()
    #     tmpi = np.array(np.arange(len(tmp)).tolist())
    #     hoge = np.argmin(np.abs(fr-f[tmp]))
    #     ifr = tmpi[tmp][hoge]
    #     fr_if = f[ifr] - klist.sg_freq # IF, Hz
    #     fr_if_kHz = int(fr_if/1e3 + 0.5) # IF, kHz
    #     kldict_new['kids_freqs'][i] = fr_if_kHz/1e3 # IF, MHz

    #     swp_data[ki].fitIQ(**fitconf)
    #     fr = swp_data[ki].fitresult.params['fr'].value # RF, Hz
    #     fr_if = fr - klist.sg_freq # IF, Hz
    #     fr_if_kHz = int(fr_if/1e3 + 0.5) # IF, kHz
    #     kldict_new['kids_freqs'][i] = fr_if_kHz/1e3 # IF, MHz
            
    # # 'blinds_freqs` will be added if only 'blinds_relfreqs' is assigned in kidslist
    # is_newblinds = ('blinds_relfreqs' in kldict_new) and (not 'blinds_freqs' in kldict_new)
    # if is_newblinds:
    #     kldict_new['blinds_freqs'] = [0]*len(klist.kids_index)
    #     for i, ki in enumerate(klist.kids_index):
    #         kldict_new['blinds_freqs'][i] = kldict_new['kids_freqs'][i] + kldict_new['blinds_relfreqs'][i] # IF, MHz

    # swp_data = readfile_swp('rhea', 
    #                         filename=str(swp_fname),
    #                         index=-1, # all
    #                         lo=klist.sg_freq)
    
    # # 'blinds_freqs` will be added if only 'blinds_relfreqs' is assigned in kidslist
    # is_newblinds = ('blinds_relfreqs' in kldict_new) and (not 'blinds_freqs' in kldict_new)
    # if is_newblinds:
    #     kldict_new['blinds_freqs'] = [0]*len(klist.kids_index)

    # for i, ki in enumerate(klist.kids_index):
    #     swp_data[ki].fitIQ(nfwhm=3)
    #     fr = swp_data[ki].fitresult.params['fr'].value # RF, Hz
    #     fr_if = fr - klist.sg_freq # IF, Hz
    #     fr_if_kHz = int(fr_if/1e3 + 0.5) # IF, kHz
    #     kldict_new['kids_freqs'][i] = fr_if_kHz/1e3 # IF, MHz
    #     if is_newblinds:
    #         kldict_new['blinds_freqs'][i] = fr_if_kHz/1e3 + kldict_new['blinds_relfreqs'][i] # IF, MHz

    with open(klist_fname, 'w') as f:
        json.dump(kldict_new, f, indent=4)

    ############################################################################# TOD
    klist = KidsList(klist_fname) # Reload KIDs list

    rate_kHz_list = [1,100,1000] # kHz
    length_sec_list = [length/1000,1,0.1] # sec
    for rate_kHz, length_sec in zip(rate_kHz_list,length_sec_list):
        length = int(rate_kHz * length_sec * 1000)
        tod_fname = path_dir.joinpath(f'tod_{label}_{rate_kHz}KSPS.rawdata')
        tod_dt = datetime.datetime.now(tz=timezone.utc)
        if not nodb:
            ## Registration to the database
            tod = RheaTod(path=tod_fname, 
                          mtime=tod_dt,
                          ntone=len(klist.tone_freqs),
                          length=length,
                          rate=rate_kHz*1e3,
                          duration=timedelta(seconds=length_sec),
                          sg=klist.sg_freq)
            tod.save()
            meas.tods.add(tod)
            meas.save()

        # todinfo = TodInfo(max_ch=fpga.MAX_CH,
        #                   dds_f_MHz=klist.tone_freqs,
        #                   data_length=length,
        #                   rate_kSPS=rate_kHz,
        #                   start_dt=tod_dt,
        #                   power=power,
        #                   amps=None,
        #                   phases=None)

        measure_tod(fpga        = fpga,
#                    MAX_CH      = MAX_CH, # for rhea_comm_old_202209
                    max_ch      = max_ch,            
#                    dds_f_MHz   = klist.tone_freqs, # for rhea_comm_old_202209
                    dds_f_megahz   = klist.tone_freqs,
                    data_length = length,
#                    rate_kSPS   = rate_kHz, # for rhea_comm_old_202209
                    rate_ksps   = rate_kHz,                    
                    power       = power,
                    fname       = tod_fname,
                    swap_dac    = swap_dac,
                    swap_adc    = swap_adc,
                    amps        = klist.tone_amps,
                    phases      = klist.tone_phases)

        # measure_tod(fpga    = fpga,
        #             todinfo = todinfo,
        #             fpath   = Path(tod_fname),
        #             verbose = True)


if __name__=='__main__':
    import os
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="config file", type=str)
    parser.add_argument("length", help="data length [sec] at 1kSPS (note: 1sec with 100kSPS, 0.1sec with 1000kSPS)", type=int)    
    parser.add_argument("-c", "--comment", help="comment", type=str, default='')
    parser.add_argument("-p", "--path_dir", help="path to the directory to store data", type=str, default=None)
    parser.add_argument("--nodb", help="Not adding database", action='store_true', default=False)    
    parser.add_argument("--section", help="section name to be run", type=str, default=None)    
    args = parser.parse_args()

    ini = make_config(args.config)

    path_dir = args.path_dir
    if path_dir is None:
        path_dir = path_creator(DATA_DIR)
    else:
        path_dir = Path(path_dir)
        path_checker(path_dir)

    from multiprocessing import Process
    import time
    for isect in ini.keys():
        if args.section is not None:
            if isect != args.section:
                print(f"{isect} is skipped.")
                continue
        confitems = ini[isect]
        print(f'FPGA: {confitems["ip_address"]}')
        try:
#            confitems['fpga'] = fpga_control(ip_address=confitems['ip_address']) # for rhea_comm_old_202209
            confitems['fpga'] = FPGAControl(ip_address=confitems['ip_address'])            
        except TimeoutError:
            print('Connection to FPGA failed.')
            print(confitems['ip_address'], 'may be invalid')
            exit(1)

        confitems['length']   = args.length * confitems['rate_kHz']*1000
        confitems['comment']  = args.comment + ' ' + confitems['comment']
        confitems['path_dir'] = path_dir
        confitems['nodb']     = args.nodb
        if confitems['label'] == "": confitems['label'] = isect
        if "sg_serial" not in confitems: confitems["sg_serial"] = None # backward compatibility
        sg_config = {}
        if confitems['sg_serial'] is not None:
            sg_config['serialnum'] = int(confitems['sg_serial'])
        else:
            sg_config['channel'] = int(confitems['sg_ch'])
        # Order is important!! Take care!! <-- should be modified
        main(confitems['klist_path'],  confitems['length'],     confitems['fpga'],
             confitems['swap_dac'], confitems['swap_adc'],
             confitems['width_sweep'], confitems['step_sweep'], confitems['path_dir'],
             sg_config,                confitems['rate_kHz'],   confitems['comment'],
             confitems['nodb'],        confitems['label'])
        time.sleep(1)
