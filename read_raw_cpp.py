import numpy as np
import pickle
from datetime import datetime, timezone
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
import os
import sys

sys.path.append('/home/gb/program/analysis/gbproc_cpp/build')
import gbirdproc as gbp
sys.path.append('./mkid_pylibs/')
sys.path.append('./analyzer_db/')
# sys.path.append('./rhea_comm/')
import mkid_pylibs as klib
from rhea_comm.lib_read_rhea import *
from kidslist import KidsList

CHIPS = np.array(['3A', '2A', '3B', '1A', '2B', '1B', '220'])
DAQS = np.array(['GB01', 'GB02', 'GB03', 'GB04', 'GB05', 'GB06', 'GB07'])
SAVEDIR = './'

class read_rawdata_cpp():
    """Class for reading KID's data based on gbproc and mkid_pylibs."""
    
    def __init__(self, meas_id, log=True, saveraw=False):
        self.meas_id = meas_id
        db = gbp.MeasDB(meas_id)
        self.daq = db.swppath.split('_')[-1][:4]
        self.chip = CHIPS[np.where(DAQS == self.daq)[0][0]]
        print(f"DAQ: {self.daq}, Chip: {self.chip}")
        
        # Read data with cpp
        try:
            swp = gbp.RheaSwpReader(db.swppath, klist.lofreq)
            tod = gbp.RheaTodReader(db.todpath, klist.lofreq)
        except:
            klist = KidsList(db.klistpath)
            swp = gbp.RheaSwpReader(db.swppath, klist.sg_freq)
            tod = gbp.RheaTodReader(db.todpath, klist.sg_freq)
            klist.blinds_freq = np.array(klist.blinds_freqs)
            klist.kids_freq = np.array(klist.kids_freqs)
        
        # Create data objects
        swpset = [klib.Swpdata(ifreq, 'I-Q', (iiq.real, iiq.imag)) 
                  for ifreq, iiq in zip(swp.freq, swp.iq)]
        todset = [[klib.TODdata(tod.time, 'I-Q', (iiq.real, iiq.imag), ifreq, 
                               info={'n_rot': tod.syncnum, 'sync_off': tod.syncoff})] 
                  for ifreq, iiq in zip(tod.freq, tod.iq)]
        
        # Quality check for blind frequencies
        good_bfreq = []
        for ibind, ibfreq in zip(klist.blinds_index, klist.blinds_freq):
            if np.abs(ibfreq) > 0.1:  # Skip very low frequencies
                if not any(np.abs(np.diff(todset[ibind][0].phase)) > 6):
                    good_bfreq.append(ibfreq)
        good_bfreq = np.array(good_bfreq)
        
        # Set nearest blind tone
        binds = []
        for ifreq in klist.kids_freq:
            ind = np.where(np.min(np.abs(good_bfreq - ifreq)) == 
                          np.abs(klist.blinds_freq - ifreq))[0][0]
            binds.append(klist.blinds_index[ind])
        self.bind = binds
        
        # Create KID analysis objects
        kr = [klib.kidana_psd.PSDAnalyzer(swp=swpset[i], tod=todset[i], ctod=todset[ibind])
              for i, ibind in enumerate(binds)]
        
        # Fitting
        self._perform_fitting(kr)
        
        # Read azimuth and elevation data
        print("Reading AzEl data")
        az = gbp.get_syncaz_rhea(tod, 41, False, False)
        el = gbp.get_syncel_rhea(tod, 0, False, False)
        
        # Process log data if requested
        if log:
            print("Reading Log data")
            log = gbp.LogContainer(az.time)
            good = self._process_log_data(log, kr)
        else:
            good = np.ones(len(az.time), dtype=bool)
            good[:1000] = False  # Skip first 1000 points
            log = None
        
        # Store results
        self.kr = kr
        self.klist = klist
        self.az = az
        self.el = el
        self.rpm = az.speed
        self.time = np.array([datetime.fromtimestamp(itime, timezone.utc) for itime in az.time])
        self.good = good
        self.log = log
        self.param = get_param(self.kr)
        
        if saveraw:
            self.save_rawdata_all()
    
    def _perform_fitting(self, kr):
        """Perform fitting for all KIDs."""
        nonphase = [None] * len(kr)
        config = get_fitconfig(self.chip)
        rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, depval, skipind, guessskip = config
        
        for i, ikr in enumerate(kr):
            # Set fitting parameters
            fitini = fitinit[np.where(np.array(initind) == i)[0][0]] if i in initind else None
            rangeini = freqranges[np.where(np.array(rangeind) == i)[0][0]] if i in rangeind else [None, None]
            twokidini = twokidfitter[np.where(np.array(twokidind) == i)[0][0]] if i in twokidind else 'gaolinbg'
            dep = depval[np.where(np.array(depind) == i)[0][0]] if i in depind else 3
            skip = guessskip[np.where(np.array(skipind) == i)[0][0]] if i in skipind else False
            
            print(f'======== Fit KID{i:02} ==========')
            ikr.fitIQ(nfwhm=-1, frqrange=rangeini, fitter=twokidini, init=fitini, dep=dep, guess_skip=skip)
            nonphase[i] = 2 * np.tan(ikr.tod.rwmdata.corphase / 2.)
        
        self.phase = nonphase
    
    def _process_log_data(self, log, kr):
        """Process log data to determine good data points."""
        goodinds = log.goodIndex(0.35, 80)  # Dome, detector temp, humidity criteria
        good = np.array([False] * (goodinds[-1] + 1))
        isok = True
        for i, j in zip(goodinds[:-1], goodinds[1:]):
            good[i:j+1] = isok
            isok = not isok
        
        # Remove glitches
        for i, ikr in enumerate(kr):
            thre = 1
            glitch_ind = np.where(np.abs(np.ediff1d(ikr.tod.rwmdata.corphase, to_end=0)) > thre)[0]
            for iglitch_ind in glitch_ind:
                if iglitch_ind < len(good):
                    good[iglitch_ind] = False
        
        return good
    
    def save_rawdata_all(self):
        """Save all raw data to pickle file."""
        ret = {
            'swp_param': self.param.to_dict(),
            'utime': self.az.time[self.good],
            'el': self.el.angle[self.good],
            'az': self.az.angle[self.good],
            'phase': {f'kid{i:02}': self.phase[i][self.good] for i in range(len(self.kr))}
        }
        
        if not os.path.isdir(SAVEDIR + 'raw_data'):
            os.mkdir(SAVEDIR + 'raw_data')
        
        filename = f'{SAVEDIR}raw_data/{self.chip}_{self.meas_id}.pkl'
        with open(filename, 'wb') as f:
            pickle.dump(ret, f)
        print(f"Raw data saved to {filename}")

def get_param(kr):
    """Extract parameters from KID analysis results."""
    df = pd.DataFrame()
    params = [[], [], [], []]  # fr, Qr, Qc, Qi
    params_str = ['fr', 'Qr', 'Qc', 'Qi']
    
    for idata in kr:
        fitparams = idata.swp.fitresult.fitparamdict
        if 'fr' in fitparams:
            for iparam, istr in zip(params, params_str):
                iparam.append(fitparams[istr])
        elif 'Qi1' in fitparams:
            for iparam, istr in zip(params, params_str):
                iparam.append(fitparams[istr + '1'])
        else:
            for iparam, istr in zip(params, params_str):
                iparam.append(fitparams[istr + '2'])
    
    for iparam, istr in zip(params, params_str):
        df[istr] = iparam
    return df


def get_fitconfig(chip):
    # configuration for fitting each detectors
    # See mkid_pylibs_sueno for each configuration.
    if True:
        if chip == '1A':
            rangeind = [15]
            freqranges = [(4.886e9, 4.8885e9)]
            twokidind = [5, 6]
            twokidfitter = ['gaolinbg2f', 'gaolinbg2l']
            initind = []
            fitinit = []
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '1B':
            rangeind = [10, 17, 19, 20, 21, 22] # KID2,3,4 are degenerate.
            freqranges = [(4.8105e9, 4.8125e9), (4.878e9, 4.881e9), (4.888e9, 4.890e9), (4.890e9, 4.892e9), (4.894e9, 4.8964e9), (4.8965e9, 4.899e9)]
            twokidind = [9, 17, 18]
            twokidfitter = ['gaolinbg2f', 'gaolinbg2f', 'gaolinbg2l']
            initind = []
            fitinit = []
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '2A':
            rangeind = []
            freqranges = []
            twokidind = []
            twokidfitter = []
            initind = []
            fitinit = []
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '2B':
            rangeind = []
            freqranges = []
            twokidind = []
            twokidfitter = []
            initind = []
            fitinit = []
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '3A':
            rangeind = [7, 8, 9, 10, 13, 14]
            freqranges = [(4.736e9, 4.7385e9), (4.736e9, 4.7385e9), (4.7385e9, 4.741e9), (4.7385e9, 4.741e9), (4.7825e9, 4.785e9), (4.7825e9, 4.785e9)]
            #twokidind = [7, 8, 9, 10, 13, 14]
            #twokidfitter = ['gaolinbg2f', 'gaolinbg2l', 'gaolinbg2f', 'gaolinbg2l', 'gaolinbg2f', 'gaolinbg2l']
            twokidind = [7, 8, 10, 13, 14]
            twokidfitter = ['gaolinbg2f', 'gaolinbg2l', 'gaolinbg2l', 'gaolinbg2f', 'gaolinbg2l']
            #initind = [7, 8, 9, 10, 13, 14]
            initind = [7, 8, 10, 13, 14]
            fitinit = [{'fr1' : 4.7373e9, 'fr2' : 4.7376e9, 'Qr1': 15000, 'Qr2': 30000, 'Qc1': 15000, 'Qc2': 30000, 'arga': 0, 'absa': 0.01, 'tau': 0, 'phi01': 0, 'phi02': 0, 'c': 0 }, # 7
                       {'fr1' : 4.7373e9, 'fr2' : 4.7376e9, 'Qr1': 15000, 'Qr2': 30000, 'Qc1': 15000, 'Qc2': 30000, 'arga': 0, 'absa': 0.01, 'tau': 0, 'phi01': 0, 'phi02': 0, 'c': 0 }, # 8
                       #{'fr1' : 4.7394e9, 'fr2' : 4.7399e9, 'Qr1': 10000, 'Qr2': 30000, 'Qc1': 17000, 'Qc2': 40000, 'arga': 0, 'absa': 0.01, 'tau': 0, 'phi01': 0, 'phi02': 0, 'c': 0 }, # 9
                       {'fr1' : 4.7395e9, 'fr2' : 4.7399e9,'Qr1': 16000, 'Qr2': 40000, 'Qc1': 17000, 'Qc2': 40000}, # 10
                       {'fr1' : 4.7836e9, 'fr2':4.7842e9, 'Qr1': 15000, 'Qr2': 35000, 'Qc1': 10000, 'Qc2': 20000}, # 13
                       {'fr1' : 4.7836e9, 'fr2':4.7842e9, 'Qr1': 15000, 'Qr2': 35000, 'Qc1': 10000, 'Qc2': 20000}] # 14
            depind = [8, 9]
            dep = [2, 2]
            #skipind = [7, 8, 9]
            #guessskip = [True, True, True]
            skipind = [7]
            guessskip = [True]
            
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '3B':
            #rangeind = [0, 12]
            #freqranges = [(4.6858e9, 4.6868e9), (4.819e9, 4.821e9)]
            rangeind = [0]
            freqranges = [(4.68575e9, 4.6868e9)]

            twokidind = []
            twokidfitter = []
            initind = [0]
            fitinit = [{'fr' : 4.6862e9, 'Qr' : 15000, 'Qc' : 6000, 'phi0' : 0.1}]
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '220':
            rangeind = [0, 2, 15] 
            freqranges = [(4.854e9, 4.857e9), (4.867e9, 4.869e9), (4.993e9, 4.995e9)]
            twokidind = [14]
            twokidfitter = ['gaolinbg2f']
            initind = []
            fitinit = []
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        else:
            print('Error : chip:' + chip + 'is not matched.')