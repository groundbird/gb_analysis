#!/usr/bin/env python3

from os  import listdir
from os.path import isdir, basename
from math import pi
from .rotary_reader import rot_reader, find_latest_file
from datetime import datetime
import numpy as np
from re import sub
import glob



datadir_pi = '/home/pi/logger/data/'
datadir_gb = '/home/gb/logger/data/'
if isdir(datadir_pi):
    datadir_head = datadir_pi
elif isdir(datadir_gb):
    datadir_head = datadir_gb

datadir_type = {}
datadir_type['thermo'] = 'thermo'
datadir_type['az'] = 'rotation_encoder'
datadir_type['accl'] = 'accl'

datatype_fileformat = {}
datatype_fileformat['thermo'] = '%Y/%m/%Y%m%d_detector.cal'
datatype_fileformat['az']     = '%Y/%m/%d/rot_%Y-%m%d-%H%M%S.dat.xz'
datatype_fileformat['accl']   = '%Y/%m/%Y%m%d.dat'

def read_date(st):
    fmt = '%Y%m%d%H%M%S.%f'
    dst = st.strip()
    now = datetime.now()
    if len(dst.split('.')[0]) == 4: dst = dst[0:4] + '00' + dst[4:]
    if len(dst.split('.')[0]) == 6: dst = now.strftime('%m%d') + dst
    if len(dst.split('.')[0]) == 10: dst = now.strftime('%y') + dst
    if len(dst.split('.')[0]) == 12: dst = now.strftime('%Y')[0:2] + dst
    if len(dst.split('.')) == 1: dst = dst.split('.')[0] + '.000'

    return datetime.strptime(dst, fmt)

def select_region(data_list, key1, key2):
    _data_list = list(data_list) + [key1, key2]
    _data_list.sort()
    num1 = _data_list.index(key1)
    num2 = _data_list.index(key2)
    cnt2 = _data_list.count(key2)
    return num1, num2 + cnt2 - 2

def conv_time(t_offset, t):
    return float(t_offset + t / 1000.)

def conv_angle(v): # raw to [rad]
    return float(2. * pi * v / 8192.)

def print_line(t, v):
    tt = t
    vv = v
    if type(tt) == float: tt = '%.4f'  % tt
    if type(vv) == float: vv = '%8.6f' % vv
    print(tt, vv)
    return

def output_file(ftype, fname,
                output_return = True,
                output_stdout = False,
                time_raw = False,
                val_raw = False):
    ret = None
    if output_return: ret = []
    if ftype == 'az':
        f = rot_reader(fname)
        t_offset = float(f.header['unixtime'])
    elif ftype == 'thermo':
        f = np.genfromtxt(fname,unpack=True)
        f = zip(f[1],f[2])
        t_offset = 0
    elif ftype == 'accl':
        f = np.genfromtxt(fname,unpack=True)
        f = zip(f[1],[f[2],f[3],f[4]])
        t_offset = 0

    for t, v in f:
        tt = t
        vv = v
        if not time_raw:
            if ftype == 'az': tt = conv_time(t_offset, t)
            elif ftype == 'thermo'  : tt = t + t_offset
            elif ftype == 'accl': tt = t + t_offset
        if not val_raw:
            if ftype == 'az': vv = conv_angle(v)
            elif ftype == 'thermo'  : vv = v
            elif ftype == 'accl': vv = v
        if output_return: ret += [(tt, vv)]
        if output_stdout: print_line(tt, vv)
        pass

    return ret

def output_date(ftype, start_st, stop_st,
                output_return = True,
                output_stdout = False,
                time_raw = False,
                val_raw = False):
    ret = None
    if output_return: ret = []
    start = read_date(start_st)
    if output_stdout : print(f'start: {start}')
    stop  = read_date(stop_st)
    if output_stdout : print(f' stop: {stop}')

    datadir = datadir_head + datadir_type[ftype]

    dir_year = listdir(datadir)
    dir_year.sort()
    a, b = select_region(dir_year, start.strftime('%Y'), stop.strftime('%Y'))
    a = max(a-1, 0)
    dir_year = dir_year[a:b]

    dir_month = []
    for d in dir_year:
        dir_month += [d+'/'+dd for dd in listdir(datadir + '/' + d)]
        pass
    dir_month.sort()
    a, b = select_region(dir_month, start.strftime('%Y/%m'), stop.strftime('%Y/%m'))
    a = max(a-1, 0)
    dir_month = dir_month[a:b]

    if datatype_fileformat[ftype].startswith('%Y/%m/%d/'):
        dir_day = []
        for d in dir_month:
            dir_day += [d+'/'+dd for dd in listdir(datadir + '/' + d)]
            pass
        dir_day.sort()
        a, b = select_region(dir_day, start.strftime('%Y/%m/%d'), stop.strftime('%Y/%m/%d'))
        a = max(a-1, 0)
        dir_day = dir_day[a:b]
        dirdates = dir_day
    else:
        dirdates = dir_month

    fname_format = datatype_fileformat[ftype]
    fname_key = basename(datatype_fileformat[ftype])
    fname_key = sub(r'%.', "*", fname_key)
    fname_key = sub(r'\*+', "*", fname_key)
    fnames = []
    for d in dirdates:
        fnames += glob.glob(datadir + '/' + d + '/' +fname_key)
        pass
    fnames.sort()

    a, b = select_region(fnames, datadir + '/' + start.strftime(fname_format), datadir + '/' + stop.strftime(fname_format))
    a = max(a-1, 0)
    fnames = fnames[a:b]

    t_start = float(start.strftime('%s.%f'))
    t_stop  = float(stop .strftime('%s.%f'))

    for fname in fnames:
        print(fname)
        if ftype == 'az':
            f = rot_reader(fname)
            t_offset = float(f.header['unixtime'])
            skip = int(round((t_start - t_offset) * 1000 - f.read(1)[0][0] - 0.5))
            if skip > 0: f.seek(skip)
        elif ftype == 'thermo':
            f = np.genfromtxt(fname,unpack=True)
            t_offset = 0
            skip = int(round((t_start - t_offset) - f[1][0] - 0.5)/10)
            #skip = max(skip,1)
            skip = 1
            f = zip(f[1][skip-1:],f[2][skip-1:])
        elif ftype == 'accl':
            f = np.genfromtxt(fname,unpack=True)
            f = zip(f[1],[f[2],f[3],f[4]])
            t_offset = 0
            skip = int(round((t_start - t_offset) - f[1][0] - 0.5)/10)
            #skip = max(skip,1)
            skip = 1
            f = zip(f[1][skip-1:],[f[2][skip-1:],f[3][skip-1:],f[4][skip-1:]])

        print(t_offset, t_start)
        print(skip)
        for t, v in f:
            tt = t
            vv = v
            if not time_raw:
                if ftype == 'az': tt = conv_time(t_offset, t)
                elif ftype == 'thermo'  : tt = t + t_offset
                elif ftype == 'accl': tt = t + t_offset
            if not val_raw:
                if ftype == 'az': vv = conv_angle(v)
                elif ftype == 'thermo'  : vv = v
                elif ftype == 'accl': vv = v
            if output_return: ret += [(tt, vv)]
            if output_stdout: print_line(tt, vv)
            if time_raw:
                if ftype == 'az': tt = conv_time(t_offset, t)
                elif ftype == 'thermo'  : tt = t + t_offset
                elif ftype == 'accl': tt = t + t_offset
            if tt > t_stop: break
            pass
        pass
    return ret

