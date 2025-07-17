#!/usr/bin/env python3

from enum import Enum
from struct import unpack
from os.path import getsize
from numpy import median
import numpy as np

class PacketError(Exception):
    def __init__(self, key):
        super().__init__(key)
        self.key = key
    def print(self):
        print(f'packet error. ({self.key})')

class PacketType(Enum):
    DATA  = 0
    SYNC  = 1
    SGSWP = 2

class ReadSgSwpFile():
    def __init__(self, _path, average=True):
        self._path    = _path
        self.n_rotate = -1
        self.sync_off = 0
        self.swp_cnt  = 0
        self.buff = None
        self._i = 0
        self.packet_size = self.get_packet_size()
        self.length = int(getsize(self._path)/self.packet_size)
        self.ch_num = int((self.packet_size - 7)/7/2)
        self._f = open(self._path, 'rb')
        self.ret = [dict() for i in range(self.ch_num)]
        self.f_start, self.f_stop, self.f_step = self.read_header()
        self.sg_freq = self.f_start + self.f_step
        self.data = self.read_file()
        if average is True:
            self.f_avr = None
            self.data_avr = self.IQ_averaging()
            pass

    def __del__(self):
        if not self._f.closed:
            self._f.close()

    def __iter__(self):
        return self

    def __next__(self):
        if self._i == self.length:
            raise StopIteration()
        self.buff = self._f.read(self.packet_size)
        timestamp, data, ptype = self.read_packet()
        self._i += 1
        if ptype == PacketType.DATA:
            return timestamp, data
        elif ptype == PacketType.SYNC:
            self.n_rotate = timestamp
            self.sync_off = data[0]
            return self.__next__()
        elif ptype == PacketType.SGSWP:
            self.sg_freq = self.f_start + self.swp_cnt * self.f_step
            if self.sg_freq == self.f_stop:
                self.swp_cnt = 0
            else:
                self.swp_cnt += 1
            return timestamp, data
        pass

    def get_packet_size(self):
        ret = []
        offset = 0
        _f = open(self._path, 'rb')
        buff = _f.read(2**12)
        _f.close()
        try:
            if buff[0] in  [0xaa, 0xff]: ret += [0]
            else: raise PacketError('get-packet-size : header')
        except PacketError as _e:
            _e.print()
            pass
        while True:
            _p = buff.find(0xee) ## 0xee = footer
            if _p == len(buff) - 1: break
            if _p == -1: break
            buff = buff[_p+1:]
            offset += _p+1
            if buff[0] in [0xff, 0xf5, 0xaa]: ret += [offset]
            pass
        ret_diff = []
        _p = ret[0]
        for _n in ret[1:]:
            ret_diff += [_n - _p]
            _p = _n
            pass
        return int(median(ret_diff))

    def read_packet(self):
        buff = self.buff
        try:
            if buff[0] == 0xff:
                ptype = PacketType.DATA
            elif buff[0] == 0xf5:
                ptype = PacketType.SYNC
            elif buff[0] == 0xaa:
                ptype = PacketType.SGSWP
            else: raise PacketError('header')
            if buff[-1] != 0xee:
                raise PacketError('footer')
            time = (unpack('b', buff[1:2])[0] << (8 * 4)) + unpack('>I', buff[2:6])[0]
            data = []
            for _dn in range(2*self.ch_num):
                d_1 = unpack('b',  buff[6 + 7 * _dn : 6 + 7 * _dn + 1])[0]
                d_2 = unpack('>H', buff[6 + 7 * _dn + 1 : 6 + 7 * _dn + 3])[0]
                d_3 = unpack('>I', buff[6 + 7 * _dn + 3 : 6 + 7 * _dn + 7])[0]
                data += [(d_1 << (8 * 6)) + (d_2 << (8 * 4)) + d_3]
                pass
        except PacketError as _e:
            _e.print()
        return time, data, ptype

    def read_header(self):
        buff = self._f.read(self.packet_size)
        rate = (unpack('b', buff[1:2])[0] << (8 * 4)) + unpack('>I', buff[2:6])[0]
        f_start = unpack('>H', buff[6 + 7 + 1 : 6 + 7 + 3])[0] # MHz
        f_stop  = unpack('>H', buff[6 + 7 + 3 : 6 + 7 + 5])[0] # MHz
        f_step  = unpack('>H', buff[6 + 7 + 5 : 6 + 7 + 7])[0] # kHz
        f_step  = float(f_step / 1000.) # kHz -> MHz
        for i in range(self.ch_num):
            f_1 = unpack('b', buff[6 + 7 * 2*i     : 6 + 7 * 2*i + 1])[0]
            f_2 = unpack('>H', buff[6 + 7 * 2*i + 1 : 6 + 7 * 2*i + 3])[0]
            f_3 = unpack('>I', buff[6 + 7 * 2*i + 3 : 6 + 7 * 2*i + 7])[0]
            freq = (f_1 << (8 * 6)) + (f_2 << (8 * 4)) + f_3
            self.ret[i]['name'    ] = self._path.split('/')[-1]
            self.ret[i]['rate'    ] = float(rate)
            self.ret[i]['freq'    ] = float(freq)
            self.ret[i]['f_start' ] = int(f_start)
            self.ret[i]['f_stop'  ] = int(f_stop)
            self.ret[i]['f_step'  ] = f_step
            self.ret[i]['n_swp_run'] = None
            self.ret[i]['time'    ] = []
            self.ret[i]['n_rot'   ] = []
            self.ret[i]['sync_off'] = []
            self.ret[i]['sg_freq' ] = []
            self.ret[i]['I'       ] = []
            self.ret[i]['Q'       ] = []
            pass
        self._i += 1
        return f_start, f_stop, f_step

    def read_file(self):
        for _ts, _d in self:
            for i in range(self.ch_num):
                self.ret[i]['time'    ].append(_ts)
                self.ret[i]['n_rot'   ].append(self.n_rotate)
                self.ret[i]['sync_off'].append(self.sync_off)
                self.ret[i]['sg_freq' ].append(self.sg_freq)
                self.ret[i]['I'       ].append(_d[2*i])
                self.ret[i]['Q'       ].append(_d[2*i+1])
                pass
            pass
        return self.ret

    def detect_step(self):
        _f = self.data[0]['sg_freq']
        _start = np.argwhere(np.diff(_f) < 0)[:, 0] + 1
        if _start.shape[0] > 1:
            _end = _start[1:] - 1
            _start = _start[:-1]
        elif _start.shape[0] == 1:
            _end = np.argwhere(np.diff(_f) != 0)[:, 0][-1] + 101
            _end = np.array([_end])
        return _start, _end, _start.shape[0]

    def f_averaging(self, _v, _start, _end, _off=0):
        _f = self.data[0]['sg_freq'][_start:_end]
        _v = _v[_start:_end]
        swp_intv = np.argwhere(np.diff(_f) != 0.)[:, 0] + 1
        swp_intv = np.insert(swp_intv, 0, 0)
        nstep = swp_intv.shape[0]
        swp_intv = np.append(swp_intv, -1)
        f_avr = []
        v_avr = []
        for j in range(nstep):
            f_avr.append(_f[swp_intv[j]+1])
            v_avr.append(np.mean(_v[swp_intv[j]+_off:swp_intv[j+1]-1]))
            pass
        if self.f_avr is None:
            self.f_avr = np.array(f_avr)
            pass
        return v_avr

    def IQ_averaging(self):
        start, end, nrun = self.detect_step()
        for i in range(self.ch_num):
            self.ret[i]['I_avr'] = []
            self.ret[i]['Q_avr'] = []
            for _r in range(nrun):
                self.ret[i]['I_avr'].append(self.f_averaging(self.data[i]['I'], start[_r], end[_r], _off=2))
                self.ret[i]['Q_avr'].append(self.f_averaging(self.data[i]['Q'], start[_r], end[_r], _off=2))
                pass
            self.ret[i]['swp_freq'] = self.f_avr + self.data[i]['freq']/1e6
            self.ret[i]['n_swp_run'] = int(nrun)
            pass
        return self.ret

