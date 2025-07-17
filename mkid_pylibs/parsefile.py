import numpy as np
import re
import datetime
import sys

def parse_riken_csv(filename):
    f, db, ang = np.loadtxt(filename, unpack=True, delimiter=',')
    return None, (f, db, ang)

def parse_vna_csv(filename):
    state = 'START'
    data = []
    info = dict()
    timefmt = '!Date: %A, %B %d, %Y %H:%M:%S'
    for i, line in enumerate(open(filename)):
        line = line.rstrip()
        if line == '':
            continue
        if line[0] == '!':
            if re.match('!Date:', line):
                info['date'] = datetime.datetime.strptime(line, timefmt)
            elif line[1:8] =='CW Freq': # read CW ferq(Hz) from csv file
                info['CW Freq'] = float(line[10:20])
                continue
        if state == 'START':
            if line[0:5] == 'BEGIN':
                state = 'HEADER'
        elif state == 'HEADER':
            info['header'] = line.split(',')
            state = 'BODY'
        elif state == 'BODY':
            if line[0:3] == 'END':
                state = 'END'
                break
            data.append([float(val) for val in line.split(',')])
        else:
            raise RuntimeError(f'ERROR:: parse error, state == {state}')
    print(info['header'])
    if tuple(info['header']) != ('Freq(Hz)', 'S21(DB)', 'S21(DEG)') and tuple(info['header']) != ('Time(s)', 'S21(DB)', 'S21(DEG)'):
        print('seems not %s data' % ('Freq(Hz)', 'S21(DB)', 'S21(DEG)'), file=sys.stderr)
    return info, np.array(data).T

def parse_vna_citi(filename):
    state = 'HEADER'
    var = []
    data = []
    info = dict()
    timefmt = '!Date: %A, %B %d, %Y %H:%M:%S'
    for i, line in enumerate(open(filename)):
        line = line.rstrip()
        if line == '':
            continue
        if line[0] == '!':
            if re.match('!Date:', line):
                info['date'] = datetime.datetime.strptime(line, timefmt)
                continue
        elif state == 'HEADER':
            if re.match('CITIFILE\s', line):
                pass
            elif re.match('NAME\s', line):
                info['NAME'] = line.split(' ')[1]
            elif re.match('VAR\s', line):
                info['VAR'] = line.split(' ')[1]
            elif re.match('DATA\s', line):
                info['DATA'] = line.split(' ')[1:]
            elif re.match('VAR_LIST_BEGIN', line):
                state = 'VARLIST'
            elif re.match('BEGIN', line):
                state = 'BODY'
        elif state == 'VARLIST':
            if line == 'VAR_LIST_END':
                state = 'HEADER'
            else:
                var.append(float(line))
        elif state == 'BODY':
            if line[0:3] == 'END':
                state = 'END'
                break
            data.append([float(val) for val in line.split(',')])
        else:
            raise RuntimeError('ERROR:: parse error, state == {state}')
    return info, np.concatenate([np.array([var]), np.array(data).T])

