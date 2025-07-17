#!/usr/bin/env python3

import subprocess
import serial
import re
from time import sleep
from argparse import ArgumentParser
import sys

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

def str_to_hex(f_str):
    result = re.search(r'\d+(\.\d+)?', f_str)    
    f_num_str = result.group()
    unit_str = f_str[result.end():]
    return frequency_formatter(float(f_num_str), getattr(FUnit, unit_str))

def hex_to_freq(f_hex, unit=FUnit.GHz):
    mhz_val = int(f_hex, 16)
    return mhz_val/unit*1e-3

def frequency_formatter(f_num, hz_scale):
    return hex_conv(f_num*hz_scale/FUnit.mHz)

def hex_conv(f_mHz, n_byte=6): # f_num should be described in mHz
    if f_mHz > 0xffffffffffff:
        raise Exception('over 12 characters')
    fmt_str = '{' + f':0{2*n_byte}X' + '}'
    return fmt_str.format(int(f_mHz))


'''
ser = serial.Serial('/dev/ttyACM0')
ser.write('0C'+str_to_hex('650MHz'))
ser.close()
id_s = '0010007FB01D1520200297'
f_test = '00AE9F7BCC00\n'
'''
class QuickSyn:
    def __init__(self, port):
        self._ser = serial.Serial(port, timeout=0.1)

    def _wr(self, command):
        self._ser.write(command.encode('utf-8'))
        sleep(0.1)
        return self._ser.readline()

    def get_id(self):
        id_str = self._wr('01')
        return ID(id_str)

    def get_status(self):
        status_str = self._wr('02')
        return Status(status_str)

    def get_frequency(self, unit=FUnit.GHz):
        readstr = self._wr('04')
        return hex_to_freq(readstr, unit=unit)

    def get_temperature(self):
        readstr = self._wr('10')
        return int(readstr, 16)/10.

    def set_freq_mHz(self, freq_mHz):
        '''sets frequency of quicksyn to given value expressed in [mHz] integer'''
        self._wr('0C'+hex_conv(freq_mHz))

    def set_freq_str(self, freq_str):
        '''sets frequency of quicksyn to given value expressed in str like 6.015GHz'''
        self._wr('0C'+str_to_hex(freq_str))

    def set_rfout(self, on=True):
        if on:
            self._wr('0F01')
        else:
            self._wr('0F00')

    def normal_sweep(self):
        f_start = '4.9GHz'
        f_stop = '5.1GHz' # span = 200MHz
        f_step = '10kHz'
        dwell_us = 1*1e5
        n_run = 10
        trigger = 0
        direction = 0

        swp_comm = ['1C',
                    str_to_hex(f_start),
                    str_to_hex(f_stop),
                    str_to_hex(f_step),
                    '0000',      # must be 0
                    hex_conv(dwell_us, n_byte=4),
                    hex_conv(n_run, n_byte=2),
                    hex_conv((1<<2)*trigger | direction, n_byte=1)]
        return self._wr("".join(swp_comm))

    def fast_sweep(self):
        f_start = '4.999GHz'
        f_stop = '5.001GHz' # span = 2MHz
        n_point = 200 # step = 10kHz
        dwell_us = 10*1e3
        n_run = 10
        trigger = 0
        direction = 2

        swp_comm = ['17',
                    str_to_hex(f_start),
                    str_to_hex(f_stop),
                    hex_conv(n_point, n_byte=2),
                    '0000',      # must be 0
                    hex_conv(dwell_us, n_byte=4),
                    hex_conv(n_run, n_byte=2),
                    hex_conv((1<<2)*trigger | direction, n_byte=1)]
        return self._wr("".join(swp_comm))

    def close(self):
        self._ser.close()


if __name__ == '__main__':
    desc = '{0} [Args] [Options]\nDetailed options -h or --help'.format(__file__)
    parser = ArgumentParser(description=desc)

    parser.add_argument('-swp', '--sweep',
                        type=str,
                        dest='sweep',
                        default=None,
                        help='Select [n]:normal_sweep or [f]:fast_sweep')

    parser.add_argument('-f', '--frequency',
                        type=str,
                        dest='frequency',
                        default=None,
                        help='Frequency with unit.\nex) 650MHz, 4.5GHz')

    parser.add_argument('--port',
                        type=str,
                        dest='port',
                        default='/dev/ttyS5',
                        help='path of COM port.')

    parser.add_argument('-p', '--power',
                        type=str,
                        dest='power',
                        default=None,
                        help='Select [on] or [off]')

    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Verbose mode.')

    args = parser.parse_args()

    # Check
    if args.power:
        if (args.power != 'on') and (args.power !='off'):
            parser.print_help()
            sys.exit()

    # Initialize QuickSyn
    quicksyn = QuickSyn(port=args.port)
    if args.verbose:
        idf = quicksyn.get_id()
        print('Model #\t\t', idf.model_number)
        print('Option #\t', idf.option_number)
        print('Soft ver\t', idf.soft_version)
        print('Serial #\t', idf.serial_number)
        print()

        status = quicksyn.get_status()
        print('External Ref\t', status.ext_ref_detected)
        print('RF locked\t', status.rf_locked)
        print('Ref locked\t', status.ref_locked)
        print('RF output\t', 'On' if status.rf_output_on else 'Off')
        print('Voltage OK\t', status.voltage_ok)
        print('Ref output\t', 'On' if status.ref_output_on else 'Off')
        print('Lock recovery\t', 'On' if status.lock_recovery else 'Off')
        print()

        freq = quicksyn.get_frequency(unit=FUnit.GHz)
        print('Current freq.\t', freq, 'GHz')
    

    if args.frequency:
        if args.verbose:
            print('Frequency set mode')
        quicksyn.set_freq_str(args.frequency)
        sleep(0.1)
        if args.verbose:
            print('Frequency set to', quicksyn.get_frequency(unit=FUnit.GHz), 'GHz')

    if args.power:
        if args.verbose:
            print('RF output will be set to', args.power)
        if args.power == 'on':
            quicksyn.set_rfout(on=True)
        elif args.power == 'off':
            quicksyn.set_rfout(on=False)
        else: # will not reach here
            pass
        if args.verbose:
            status = quicksyn.get_status()
            print('RF output was set to', 'On' if status.rf_output_on else 'Off')

    if args.sweep:
        if args.sweep == 'n':
            ret_swp = quicksyn.normal_sweep()
        elif args.sweep == 'f':
            ret_swp = quicksyn.fast_sweep()
        else:
            print('invalid argument for -swp')

    if (args.frequency is None) and (args.power is None) and (args.sweep is None):
        if not args.verbose:
            freq = quicksyn.get_frequency(unit=FUnit.GHz)
            print('Current freq.\t', freq, 'GHz')
            parser.print_help()

    quicksyn.close()
