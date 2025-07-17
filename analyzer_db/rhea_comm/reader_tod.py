#!/usr/bin/env python3

from packet_reader import read_file, get_length
from sys   import argv, stderr
from numpy import mean
from os.path import isfile

try:
    args = argv[1:]
    arg_val = []
    mode = 'normal'
    isNormalize = False
    while args:
        if args[0] == '-h':
            mode = 'header'
            args = args[1:]
        elif args[0] == '-r':
            mode = 'raw'
            args = args[1:]
        elif args[0] == '-n':
            isNormalize = True
            args = args[1:]
        else:
            arg_val += [args[0]]
            args = args[1:]
            pass
        pass
    if isNormalize and mode == 'raw':
        print('-r and -n are conflicted.', file=stderr)
        raise
    offset = 0
    length = None
    if len(arg_val) < 1: raise
    fname = arg_val[0]
    if len(arg_val) >= 2: offset += int(arg_val[1])
    if len(arg_val) >= 3: length  = max(int(arg_val[2]) - offset, 0)

    if not isfile(fname):
        print(f'{fname:s} is not existed.', file=stderr)
        raise
except Exception as e:
    print('Usage: reader_tod.py [filename] ([start] [stop])', file=stderr)
    print('       option: -h: output of measurement setting (header)', file=stderr)
    print('               -r: output of raw_data with header', file=stderr)
    print('               -n: normalize to 1-ch-readout', file=stderr)
    exit(1)
    pass


## tool
def normal_read():
    for time, data, tmp, tmp1 in read_file(fname, length = length, offset = offset, sync=True):
        print(time, end=' ')
        for d in data: print(d * norm_factor, end=' ')
        print()
        pass
    return

def header_read():
    time = 0
    data = 0
    for t, d, n_r, s_off in read_file(fname, length = 0):
        time = t
        data = [v for i, v in enumerate(d) if i % 2 == 0]
        pass
    return time , data

def header_read_snap():
    ret = []
    for t, d in read_file(fname):
        if t >= 0: break
        ret += [d[0]]
        pass
    return ret

def setting_read():
    rate_SPS, freq_Hz = header_read()
    if rate_SPS == 0:
        print( 'It is old-version file.')
        print(f'length: {get_length(fname):d})')
        exit(0)
        pass
    if rate_SPS > 0:
        print(f'rate: {float(rate_SPS)/1000:f} kSPS')
    else:
        print('It is snapshot file.')
        freq_Hz = header_read_snap()
        pass
    print(f'channel: {(len(freq_Hz)):d} ch')
    for i, v in enumerate(freq_Hz):
        print(f'ch{i:03d}: {(float(v) / 1e6):7.3f} MHz')
        pass
    print(f'length: {(get_length(fname)-1):d}')
    return


## main
if mode == 'header':
    setting_read()
    exit(0)
    pass

## 'normal', 'snap', 'raw'
norm_factor = 1
if mode != 'raw':
    rate_SPS, freq_Hz = header_read()
    if rate_SPS ==  0: mode = 'raw' # old file
    if rate_SPS == -1: mode = 'snap'
    if mode == 'normal' and isNormalize:
        norm_factor = len(freq_Hz) * (rate_SPS/1000.)
        pass
    pass

if mode == 'normal': offset += 1 ## skip header packet
if mode == 'snap'  : offset += len(header_read_snap())

normal_read()
