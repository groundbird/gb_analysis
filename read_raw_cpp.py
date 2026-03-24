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


# Import configuration from separate file
from config import CHIPS, DAQS, SAVEDIR, get_fitconfig
from auto_fitconfig import get_auto_fitconfig

import socket



class read_rawdata_cpp():
    """Class for reading KID's data based on gbproc and mkid_pylibs.
    Args:
    meas_id (int): Measurement ID to process
    log (bool): Enable environmental condition filtering
    saveraw (bool): Save processed data to pickle file
    """

    def __init__(self, meas_id, log=True, saveraw=False, auto_fitconfig=False):
        self.auto_fitconfig = auto_fitconfig
        self.meas_id = meas_id
        db = gbp.MeasDB(meas_id)
        swppath   = self._fix_path(db.swppath)
        todpath   = self._fix_path(db.todpath)
        klistpath = self._fix_path(db.klistpath)

        self.daq = swppath.split('_')[-1][:4]
        self.chip = CHIPS[np.where(DAQS == self.daq)[0][0]]
        print(f"DAQ: {self.daq}, Chip: {self.chip}")

        # Read data with cpp
        try:
            swp = gbp.RheaSwpReader(swppath, klist.lofreq)
            tod = gbp.RheaTodReader(todpath, klist.lofreq)
        except:
            klist = KidsList(klistpath)
            swp = gbp.RheaSwpReader(swppath, klist.sg_freq)
            tod = gbp.RheaTodReader(todpath, klist.sg_freq)
            klist.blinds_freq = np.array(klist.blinds_freqs)
            klist.kids_freq = np.array(klist.kids_freqs)

        # Create data objects
        swpset = [klib.Swpdata(ifreq, 'I-Q', (iiq.real, iiq.imag))
                  for ifreq, iiq in zip(swp.freq, swp.iq)]
        todset = [[klib.TODdata(tod.time, 'I-Q', (iiq.real, iiq.imag), ifreq,
                               info={'n_rot': tod.syncnum, 'sync_off': tod.syncoff})]
                  for ifreq, iiq in zip(tod.freq, tod.iq)]
        self.swpset = swpset
        self.todset = todset

        # Quality check for blind frequencies
        valid_mask = np.abs(klist.blinds_freq) > 0.1
        valid_bind  = np.array(klist.blinds_index)[valid_mask]
        valid_bfreq = klist.blinds_freq[valid_mask]

        good_bfreq = []
        good_bind  = []
        for ibind, ibfreq in zip(valid_bind, valid_bfreq):
            if not any(np.abs(np.diff(todset[ibind][0].phase)) > 6):
                good_bfreq.append(ibfreq)
                good_bind.append(ibind)
        good_bfreq = np.array(good_bfreq)
        good_bind  = np.array(good_bind)

        if len(good_bfreq) == 0:
            print("WARNING: No good blind tones found. Using all valid blind tones.")
            good_bfreq = valid_bfreq
            good_bind  = valid_bind

        # Set nearest blind tone
        binds = []
        for ifreq in klist.kids_freq:
            ind = np.argmin(np.abs(good_bfreq - ifreq))
            binds.append(good_bind[ind])
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

    def _path_changer(self,path, old='/data/gb', new='/data.kiwi/gb'):
        """Replace directory name"""
        return path.replace(old, new)

    def _fix_path(self, path):
        """Apply _path_changer only when running on moa server."""
        if socket.gethostname() == 'moa':
            return self._path_changer(path)
        return path

    def _perform_fitting(self, kr):
        """Perform fitting for all KIDs."""
        nonphase = [None] * len(kr)
        if self.auto_fitconfig:
            print("Using auto fitconfig")
            config = get_auto_fitconfig(kr)
        else:
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
            ikr._fitter = twokidini


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

 # ============ basic plot ===========
    def plot_log(self):
        fig, ax = plt.subplots(3, 1, sharex=True, figsize=(10, 8))
        ax[0].plot(self.time[self.good], self.log.thermo.sync_Tdet[self.good])
        ax[0].set_ylabel('Detector temp. [K]')
        ax[1].plot(self.time[self.good], self.log.gaulli.sync_pwv[self.good])
        ax[1].set_ylabel('PWV [mm]')
        ax[2].plot(self.time[self.good], self.log.bme.sync_humidity[self.good])
        ax[2].set_ylabel('humidity [%]')
        ax[2].set_xlabel('Date')
        fig.suptitle(f'log_{self.meas_id}_{self.daq}_{self.chip}')
        fig.tight_layout()
        dirpath = SAVEDIR + 'plot/log/'
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)
        fig.savefig(dirpath + f'log_{self.meas_id}_{self.daq}_{self.chip}.jpg')
        plt.close()

    def plot_swpamp(self):
        plt.figure(figsize=(20, 14))
        for i, ikr in enumerate(self.kr):
            plt.subplot(5, 5, i + 1)
            plt.plot(ikr.swp.f, ikr.swp.amplitude, '.:', c='b', label='data')
            try:
                s = ikr.swp.fitresult.fitrange
                plt.plot(ikr.swp.f[s], ikr.swp.fitdata.amplitude[s], '-', c='r', label='fit')
            except:
                pass
            plt.axvline(ikr.tod.f, c='green', label='TOD')

            # Plot fitted fr values
            params = ikr.swp.fitresult.params

            if 'fr' in params:
                # Single resonance: always orange
                fr_val = params['fr'].value
                idx = np.argmin(np.abs(ikr.swp.f - fr_val))
                plt.plot(fr_val, ikr.swp.amplitude[idx], '*',
                        color='orange', markersize=12, label='fr')

            elif 'fr1' in params and 'fr2' in params:
                # Double resonance: determine target fr by fitter type
                # gaolinbg2f -> fr1 (first) is the target: orange
                # gaolinbg2l -> fr2 (last)  is the target: orange
                fitter = getattr(ikr, '_fitter', 'gaolinbg')
                if fitter == 'gaolinbg2f':
                    colors = ('orange', 'cyan')
                else:  # gaolinbg2l
                    colors = ('cyan', 'orange')

                fr1_val = params['fr1'].value
                fr2_val = params['fr2'].value
                for fr_val, color, label in zip(
                        (fr1_val, fr2_val), colors, ('fr1', 'fr2')):
                    idx = np.argmin(np.abs(ikr.swp.f - fr_val))
                    plt.plot(fr_val, ikr.swp.amplitude[idx], '*',
                            color=color, markersize=12, label=label)

            plt.legend(fontsize=7)
            plt.ylabel('Amplitude [a.u.]')
            plt.xlabel('Frequency [Hz]')
            plt.title(f'kid{i:02}')
        plt.suptitle(f'swpamp_{self.meas_id}_{self.daq}_{self.chip}')
        plt.tight_layout()
        dirpath = SAVEDIR + 'plot/swpamp/'
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)
        plt.savefig(dirpath + f'swpamp_{self.meas_id}_{self.daq}_{self.chip}.jpg')
        plt.close()



    def plot_bswpamp(self):
        plt.figure(figsize=(20, 20))
        for i, bi in enumerate(self.bind):
            plt.subplot(5, 5, i + 1)
            plt.plot(self.swpset[bi].f, np.abs(self.swpset[bi].iq), '.:', c='b', label='data')
            plt.axvline(np.median(self.swpset[bi].f), c='green', label='TOD')
            plt.legend()
            plt.ylabel('Amplitude [a.u.]')
            plt.xlabel('Frequency [Hz]')
            plt.title(f'kid{bi:02}')
        plt.suptitle(f'bswpamp_{self.meas_id}_{self.daq}_{self.chip}')
        plt.tight_layout()
        dirpath = SAVEDIR + 'plot/bswpamp/'
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)
        plt.savefig(dirpath + f'bswpamp_{self.meas_id}_{self.daq}_{self.chip}.jpg')
        plt.close()

    def plot_swpiq(self):
        plt.figure(figsize=(30, 20))
        for i, ikr in enumerate(self.kr):
            plt.subplot(5, 5, i + 1)
            try:
                plt.plot(ikr.swp.rwdata.i, ikr.swp.rwdata.q, '.:', c='b', label='data')
                plt.plot(ikr.swp.fitdata.rwdata.i, ikr.swp.fitdata.rwdata.q, '-', c='r', label='fit')
                plt.plot(ikr.tod.rwdata.i[self.good], ikr.tod.rwdata.q[self.good], '.', c='green', label='TOD', markersize=10)
            except:
                pass
            plt.legend()
            plt.axis('equal')
            plt.ylabel('Q [a.u.]')
            plt.xlabel('I [a.u.]')
            plt.title(f'kid{i:02}')
        plt.suptitle(f'swpiq_{self.meas_id}_{self.daq}_{self.chip}')
        plt.tight_layout()
        dirpath = SAVEDIR + 'plot/swpiq/'
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)
        plt.savefig(dirpath + f'swpiq_{self.meas_id}_{self.daq}_{self.chip}.jpg')
        plt.close()

    def plot_psd(self):
        for ikr in self.kr:
            ikr.calcpsd(dofit=False)
        plt.figure(figsize=(30, 20))
        for i, ikr in enumerate(self.kr):
            plt.subplot(5, 5, i + 1)
            plt.semilogx(ikr.psd.f, 10 * np.log10(ikr.psd.amplitude), '.:', c='r', label='amp')
            plt.semilogx(ikr.psd.f, 10 * np.log10(ikr.psd.phase), '.:', c='b', label='phase')
            plt.legend()
            plt.ylabel('PSD [dBc/Hz]')
            plt.xlabel('Frequency [Hz]')
            plt.title(f'kid{i:02}')
        plt.suptitle(f'psd_{self.meas_id}_{self.daq}_{self.chip}')
        plt.tight_layout()
        dirpath = SAVEDIR + 'plot/psd/'
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)
        plt.savefig(dirpath + f'psd_{self.meas_id}_{self.daq}_{self.chip}.jpg')
        plt.close()

    def plot_tod(self):
        plt.figure(figsize=(30, 20))
        for i, bi in enumerate(self.bind):
            plt.subplot(5, 5, i + 1)
            plt.plot(self.todset[i][0].amplitude, '.:', c='r', label='kid')
            plt.plot(self.todset[bi][0].amplitude, '.:', c='b', label='blind')
            plt.legend()
            plt.ylabel('Amplitude [a.u.]')
            plt.xlabel('time')
            plt.title(f'kid{i:02}')
        plt.suptitle(f'tod_{self.meas_id}_{self.daq}_{self.chip}')
        plt.tight_layout()
        dirpath = SAVEDIR + 'plot/tod/'
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)
        plt.savefig(dirpath + f'tod_{self.meas_id}_{self.daq}_{self.chip}.jpg')
        plt.close()

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
