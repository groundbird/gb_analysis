#!/usr/bin/env python3

from packet_reader import read_file
from sys   import argv, stderr
from numpy import mean, angle
from os.path import isfile

try:
    fname = argv[1]
    if not isfile(fname):
        print(f'{fname:s} is not existed.', file=stderr)
        raise
except Exception as e:
    print('Usage: reader_swp.py [filename]', file=stderr)
    exit(1)


## tool
def print_data(freq, data):
    if not freq: return
    if not data: return
    mdata = [mean(d) for d in zip(*data)]
    # for i in range(len(freq)):
    #     print freq[i], mdata[2*i], mdata[2*i+1],
    #     pass
    # print
    #comp_d = mdata[0] + mdata[1] * 1j
    #print(freq[0], mdata[0], mdata[1], abs(comp_d), angle(comp_d))
    ch_number = len(freq)
    for i in range(ch_number):
        comp_d = mdata[i*2] + mdata[i*2+1] * 1j
        print(freq[i], mdata[i*2], mdata[i*2+1], abs(comp_d), angle(comp_d), end=' ')
    print()
    #print(freq[0], mdata[0], mdata[1], abs(comp_d), angle(comp_d))
    pass


## main
freq_buf = []
data_buf = [None,]

for time, data, tmp, tmp1 in read_file(fname, sync=True):
    if time == 0:
        if freq_buf and not data_buf: # error skip
            print(f'readout timing error: {freq_buf[0]:d}', file=stderr)
            continue
        print_data(freq_buf, data_buf)
        freq_buf = [d for i, d in enumerate(data) if i%2 == 0]
        data_buf = []
    else:
        data_buf += [data]
        pass
    pass
print_data(freq_buf, data_buf)


