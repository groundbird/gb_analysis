#!/usr/bin/env python3

import re
from time import sleep
from argparse import ArgumentParser
import sys
import serial

class FUnit(object):
    mHz, Hz, kHz, MHz, GHz = 1e-3, 1.0, 1e+3, 1e+6, 1e+9

query_dict = {'get_id': '01',
              'get_status': '02',
              'get_freq': '04',
              'get_temperature': '10'}

class Status(object):
    def __init__(self, status_str):
        bincode = bin(int(status_str[0:2], 16))[2:][::-1]
        self.ext_ref_detected = bincode[0] == '1'
        self.rf_locked        = bincode[1] == '0'
        self.ref_locked       = bincode[2] == '0'
        self.rf_output_on     = bincode[3] == '1'
        self.voltage_ok       = bincode[4] == '0'
        self.ref_output_on    = bincode[5] == '1'
        self.lock_recovery    = bincode[7] == '1'

class ID(object):
    def __init__(self, id_str):
        self.model_number  = id_str[0:4]
        self.option_number = id_str[4:8]
        self.soft_version  = id_str[8:12]
        self.serial_number = int(id_str[12:22])

def hex_to_freq(f_hex, unit=FUnit.GHz):
    mhz_val = int(f_hex, 16)
    return mhz_val/unit*1e-3

def hz_to_hex(f_hz):
    return hex_conv(f_hz/FUnit.mHz)

def hex_conv(f_mHz, n_byte=6): # f_num should be described in mHz
    if f_mHz > 0xffffffffffff:
        raise Exception('over 12 characters')
    fmt_str = '{' + f':0{2*n_byte}X' + '}'
    return fmt_str.format(int(f_mHz))

def fstr_to_hz(f_str):
    result = re.search(r'\d+(\.\d+)?', f_str)    
    f_num_str = result.group()
    unit_str = f_str[result.end():]
    f_ret = float(f_num_str) * getattr(FUnit, unit_str)
    return round(f_ret)

def fstr_to_hex(f_str):
    return hz_to_hex(fstr_to_hz(f_str))

class Sweep:
    def __init__(self, port, mode, start, stop, step, points, dwell, run):
        self._ser      = serial.Serial(port, timeout=0.1)
        self.mode      = mode
        self.f_start   = start
        self.f_stop    = stop
        self.f_step    = step
        self.points    = int(points)
        self.dwell     = int(dwell)
        self.run_times = int(run)

    def _wr(self, command):
        self._ser.write(command.encode('utf-8'))
        sleep(0.1)
        return self._ser.readline()

    def send_command(self):
        trigger = 0
        direction = 0

        if self.mode == 'normal':
            main = '1C'
            fineness = hz_to_hex(self.f_step)
        elif self.mode == 'fast':
            main = '17'
            fineness = hex_conv(self.points, n_byte=2)

        swp_comm = [main,
                    hz_to_hex(self.f_start),
                    hz_to_hex(self.f_stop),
                    fineness,
                    '0000',      # must be 0
                    hex_conv(self.dwell, n_byte=4),
                    hex_conv(self.run_times, n_byte=2),
                    hex_conv((1<<2)*trigger | direction, n_byte=1)]
        return self._wr("".join(swp_comm))

    def set_frequency(self, fstr):
        self._wr('0C'+fstr_to_hex(fstr))

    def get_id(self):
        id_str = self._wr('01')
        return ID(id_str)

    def get_status(self):
        status_str = self._wr('02')
        return Status(status_str)

    def get_frequency(self, unit=FUnit.GHz):
        readstr = self._wr('04')
        return hex_to_freq(readstr, unit=unit)

    def close(self):
        self._ser.close()


if __name__ == '__main__':
    desc = '{0} [Args] [Options]\nDetailed options -h or --help'.format(__file__)
    parser = ArgumentParser(description=desc)

    parser.add_argument('mode',
                        type=str,
                        default='normal',
                        help='Select [normal] or [fast].')

    parser.add_argument('--port',
                        type=str,
                        dest='port',
                        default='/dev/ttyS5',
                        help='path of COM port.\n')

    parser.add_argument('-fc', '--f_center',
                        type=str,
                        dest='f_center',
                        default='4GHz',
                        help='center frequency of sweep with unit.\n')

    parser.add_argument('-fw', '--f_width',
                        type=str,
                        dest='f_width',
                        default='2MHz',
                        help='sweeping frequency half-width with unit.')

    parser.add_argument('-fs', '--f_step',
                        type=str,
                        dest='f_step',
                        default='10kHz',
                        help='sweep step frequency with unit.\n')

    parser.add_argument('-p', '--points',
                        type=str,
                        dest='points',
                        default='200',
                        help='# of sweep points.\n')

    parser.add_argument('-d', '--dwell',
                        type=str,
                        dest='dwell',
                        default='1000',
                        help='dwell time in us.\n')

    parser.add_argument('-r', '--run_times',
                        type=str,
                        dest='run_times',
                        default='10',
                        help='# of run times.\n')

    args = parser.parse_args()

    if (args.mode != 'normal') and (args.mode != 'fast'):
        print('invalid argument. first argument must be \'normal\' or \'fast\'.')
        sys.exit()

    f_center_Hz = fstr_to_hz(args.f_center)
    f_width_Hz = fstr_to_hz(args.f_width)
    f_start_Hz = f_center_Hz - f_width_Hz
    f_stop_Hz  = f_center_Hz + f_width_Hz

    sweep = Sweep(port   = args.port,
                  mode   = args.mode,
                  start  = f_start_Hz,
                  stop   = f_stop_Hz,
                  step   = fstr_to_hz(args.f_step),
                  points = args.points,
                  dwell  = args.dwell,
                  run    = args.run_times)

    ret_swp = sweep.send_command()

    sweep.close()
