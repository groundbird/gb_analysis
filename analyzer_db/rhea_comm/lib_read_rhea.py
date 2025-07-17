#!/usr/bin/env python3

from .packet_reader import read_file
from numpy import mean, angle, array, log10
from enum import Enum
from datetime import datetime
from pathlib import Path

from .ReadSgSwp import ReadSgSwpFile

import numpy as np
import warnings

def insert_data(ret, freq, data, ind=0):
    if not freq: return ret
    if not data: return ret
    mdata = [mean(d) for d in zip(*data)]
    ret['freq'].append(freq[ind])
    ret['I'   ].append(mdata[ind*2])
    ret['Q'   ].append(mdata[ind*2+1])
    return ret

def read_rhea_mulswp(fname,miniret=False):
    return read_rhea_swp(fname,ismult=True,miniret=miniret)

def read_rhea_swp(fname,ismult=False,miniret=False):
    freq_buf = []
    data_buf = [None,]
    ret = None
    for time, data, tmp, tmp1 in read_file(fname, sync=True):
        if ret is None:
            ch_num = int(len(data)/2)
            ret = [dict() for i in range(ch_num)]
            for i in range(ch_num):
                ret[i]['name'] = fname
                ret[i]['freq'] = []
                ret[i]['I'   ] = []
                ret[i]['Q'   ] = []
                pass
        if time==0:
            if freq_buf and not data_buf:
                continue
            for i in range(ch_num):
                ret[i] = insert_data(ret[i], freq_buf, data_buf, ind=i)
            freq_buf = [d for i, d in enumerate(data) if i%2 == 0]
            data_buf = []
        else:
            data_buf += [data]
            pass
        pass

    for i in range(ch_num):
        ret[i] = insert_data(ret[i], freq_buf, data_buf, ind=i)
        ret[i]['freq'] = array(ret[i]['freq']).astype(float)
        ret[i]['I']    = array(ret[i]['I']).astype(float) / 200000. / (2**28)
        ret[i]['Q']    = array(ret[i]['Q']).astype(float) / 200000. / (2**28)
        if not miniret:
            ret[i]['IQ']   = ret[i]['I'] + ret[i]['Q'] * 1j
            ret[i]['amp_rad'] = abs(ret[i]['IQ'])
            ret[i]['phase'  ] = angle(ret[i]['IQ'])
            meanIQ = mean(ret[i]['IQ'])
            ret[i]['pha_rad'] = angle(ret[i]['IQ'] / meanIQ) * abs(meanIQ)
            ret[i]['ampDB' ]  = log10(ret[i]['amp_rad']) * 20
        pass

    if ismult is True:
        return ret
    else:
        return ret[0]

def read_rhea_tod(fname,nmax=None, miniret=False):
    initialize = True
    for it,(time, data) in enumerate(read_file(fname)):
        if initialize:
            initialize = False
            ch_num = int(len(data)/2)
            ret = [dict() for i in range(ch_num)]
            for i in range(ch_num):
                ret[i]['name'] = fname
                ret[i]['rate'] = float(time)
                ret[i]['freq'] = float(data[2*i])
                ret[i]['time'] = []
                ret[i]['I']    = []
                ret[i]['Q']    = []
                pass
            pass
        else:
            for i in range(ch_num):
                ret[i]['time'].append(time)
                ret[i]['I'   ].append(data[2*i])
                ret[i]['Q'   ].append(data[2*i+1])
                pass
            pass
        pass
        if nmax is not None and it>nmax: break
    for i in range(ch_num):
        ret[i]['time'] = array(ret[i]['time']).astype(float) / ret[i]['rate']
        ret[i]['I']    = array(ret[i]['I']).astype(float)
        ret[i]['I']   /= (2**28) * 200.e6 / ret[i]['rate']
        ret[i]['Q']    = array(ret[i]['Q']).astype(float)
        ret[i]['Q']   /= (2**28) * 200.e6 / ret[i]['rate']
        if not miniret:
            ret[i]['IQ']   = ret[i]['I'] + ret[i]['Q'] * 1j
            ret[i]['amp_rad'] = abs(ret[i]['IQ'])
            ret[i]['phase'  ] = angle(ret[i]['IQ'])
            meanIQ = mean(ret[i]['IQ'])
            ret[i]['pha_rad'] = angle(ret[i]['IQ'] / meanIQ) * abs(meanIQ)
            ret[i]['ampDB' ]  = log10(ret[i]['amp_rad']) * 20
        pass

    return ret

def read_rhea_tod_sync(fname, nbegin=0, nmax=None, nstep=1, miniret=False):
    initialize = True
    for it, (time, data, n_rot, offset) in enumerate(read_file(fname, sync=True, length=nmax, offset=nbegin, step=nstep)):
        if initialize:
            initialize = False
            ch_num = int(len(data)/2)
            ret = [dict() for i in range(ch_num)]
            for i in range(ch_num):
                ret[i]['name']     = fname
                ret[i]['rate']     = float(time)
                ret[i]['freq']     = float(data[2*i])
                ret[i]['time']     = []
                ret[i]['n_rot']    = []
                ret[i]['sync_off'] = []
                ret[i]['I']        = []
                ret[i]['Q']        = []
                pass
            pass
        else:
            for i in range(ch_num):
                ret[i]['time'    ].append(time)
                ret[i]['n_rot'   ].append(n_rot)
                ret[i]['sync_off'].append(offset)
                ret[i]['I'       ].append(data[2*i])
                ret[i]['Q'       ].append(data[2*i+1])
                pass
            pass
        pass

    for i in range(ch_num):
        ret[i]['time']     = array(ret[i]['time']).astype(float) / ret[i]['rate']
        ret[i]['n_rot'   ] = array(ret[i]['n_rot']).astype(int)
        ret[i]['sync_off'] = array(ret[i]['sync_off']).astype(int)
        ret[i]['I']        = array(ret[i]['I']).astype(float)
        ret[i]['I']       /= (2**28) * 200.e6 / ret[i]['rate']
        ret[i]['Q']        = array(ret[i]['Q']).astype(float)
        ret[i]['Q']       /= (2**28) * 200.e6 / ret[i]['rate']
        if not miniret:
            ret[i]['IQ']      = ret[i]['I'] + ret[i]['Q'] * 1j
            ret[i]['amp_rad'] = abs(ret[i]['IQ'])
            ret[i]['phase'  ] = angle(ret[i]['IQ'])
            meanIQ = mean(ret[i]['IQ'])
            ret[i]['pha_rad'] = angle(ret[i]['IQ'] / meanIQ) * abs(meanIQ)
            ret[i]['ampDB' ]  = log10(ret[i]['amp_rad']) * 20
        pass

    return ret


def read_rhea_sgswp(fname):
    '''
    ret includes 
    name, rate, freq, f_start, f_stop, f_step            : len=1 
    time, n_rot, sync_off, sg_freq, I, Q, IQ, amp, phase : array of len=tod_data_length
    swp_freq, *_avr                                      : array of len=#_of_swp_step.
                                                         : *_avr is array of swp run time.
    '''
    sgswp = ReadSgSwpFile(fname)
    ch_num = sgswp.ch_num
    ret = sgswp.data_avr
    for i in range(ch_num):
        ret[i]['time']     = array(ret[i]['time']).astype(float) / ret[i]['rate']
        ret[i]['n_rot']    = array(ret[i]['n_rot']).astype(int)
        ret[i]['sync_off'] = array(ret[i]['sync_off']).astype(int)
        ret[i]['sg_freq']  = array(ret[i]['sg_freq']).astype(float)
        ret[i]['I']        = array(ret[i]['I']).astype(float)
        ret[i]['I']       /= (2**28) * 200.e6 / ret[i]['rate']
        ret[i]['Q']        = array(ret[i]['Q']).astype(float)
        ret[i]['Q']       /= (2**28) * 200.e6 / ret[i]['rate']
        ret[i]['IQ']       = ret[i]['I'] + ret[i]['Q'] * 1j
        ret[i]['amp_rad']  = abs(ret[i]['IQ'])
        ret[i]['phase'  ]  = angle(ret[i]['IQ'])
        meanIQ             = mean(ret[i]['IQ'])
        ret[i]['pha_rad']  = angle(ret[i]['IQ'] / meanIQ) * abs(meanIQ)
        ret[i]['amp_dBm']  = log10(ret[i]['amp_rad']) * 20 + 10
        ret[i]['I_avr']    = array(ret[i]['I_avr']).astype(float)
        ret[i]['I_avr']   /= (2**28) * 200.e6 / ret[i]['rate']
        ret[i]['Q_avr']    = array(ret[i]['Q_avr']).astype(float)
        ret[i]['Q_avr']   /= (2**28) * 200.e6 / ret[i]['rate']
        ret[i]['IQ_avr']   = ret[i]['I_avr'] + ret[i]['Q_avr'] * 1j
        ret[i]['amp_rad_avr']  = abs(ret[i]['IQ_avr'])
        ret[i]['phase_avr'  ]  = angle(ret[i]['IQ_avr'])
        meanIQ_avr             = mean(ret[i]['IQ_avr'])
        ret[i]['pha_rad_avr']  = angle(ret[i]['IQ_avr'] / meanIQ) * abs(meanIQ)
        ret[i]['amp_dBm_avr']  = log10(ret[i]['amp_rad_avr']) * 20 + 10
        pass

    return ret

class RheaFileType(str, Enum):
    swp = 'swp'
    tod = 'tod'
        
class RheaHeader:
    def __init__(self, path):
        self.path = Path(path)
        file_dsc = read_file(self.path)
        t, data_0 = next(file_dsc)
        file_dsc.close()
        
        self.n_ch = len(data_0)/2
        self.freq = np.array(data_0[::2])
        self.b_pack = 1 + 5 + 7*len(data_0) + 1
        info = Path(path).stat()
        n_pack = info.st_size / self.b_pack
        if not n_pack.is_integer():
            warnings.warn(f'File size and packet size do not match: {info.st_size}/{self.b_pack}')
        self.n_pack = int(n_pack)

        self.type = RheaFileType.swp if t == 0 else RheaFileType.tod
        self._t = t
        self.mtime = datetime.fromtimestamp(info.st_mtime)
    
class TodHeader(RheaHeader):
    def __init__(self, path):
        if self.type != RheaFileType.tod:
            raise Exception(f'This is not TOD: {self.path}')
        self.length = self.n_pack - 1
        self._rate = self._t

    def __new__(cls, arg):
        if isinstance(arg, RheaHeader):
            arg.__class__ = cls
            return arg
        else:
            self = super().__new__(cls)
            super(cls, self).__init__(arg)
            return self

    @property
    def sampling_rate(self):
        return self._rate
    
    @property
    def duration(self):
        return self.length/self.sampling_rate
    
class SwpHeader(RheaHeader):
    N_MEAN = 10
    def __init__(self, arg):
        if self.type != RheaFileType.swp:
            raise Exception(f'This is not swp: {self.path}')
        self.length = int(self.n_pack / (1 + self.N_MEAN))
        if self.n_pack % (1 + self.N_MEAN) != 0:
            raise Exception(f'Length do not match. {self.n_pack}/{1 + self.N_MEAN}')

        ####################### Resolution finder
        file_dsc = read_file(self.path)
        for _ in range(1+self.N_MEAN):
            next(file_dsc)
        code, data_1 = next(file_dsc)
        file_dsc.close()

        if code != 0:
            raise Exception('Cannot obtain resolution')
        df = np.array(data_1[::2] - self.freq)
        if not np.all(df == df[0]):
            raise Exception(f'Resolutions do not match: {df}')
        self.resolution = df[0]

    def __new__(cls, arg):
        if isinstance(arg, RheaHeader):
            arg.__class__ = cls
            return arg
        else:
            self = super().__new__(cls)
            super(cls, self).__init__(arg)
            return self
    
    @property
    def span(self):
        return self.resolution*self.length


if __name__ == '__main__':
    from sys import argv, stderr
    import os
    bn = os.path.basename(argv[1])
    if bn.split('_')[0] == 'swp':
        d = read_rhea_swp(argv[1])
        print(f'{d["freq"][0]:10.3e} to {d["freq"][-1]:10.3e} [Hz] -- {list(d.keys())}')
    if bn.split('_')[0] == 'swpmul':
        d = read_rhea_swp(argv[1], ismult=True)
        for dd in d: print(f'{dd["freq"][0]:10.3e} to {dd["freq"][-1]:10.3e} [Hz] -- {list(dd.keys())}')
    elif bn.split('_')[0] == 'tod':
        d = read_rhea_tod(argv[1])
        for dd in d: print(f'{dd["freq"]:10.3e} Hz, {dd["rate"]:8.0f} SPS -- {list(dd.keys())}')

    else:
        print(f'Invalid fname : {argv[0]}', file=stderr)
        print(f'(should be started by \'swp_\', \'tod_\' or \'swpmul_\')')
    exit(0)
