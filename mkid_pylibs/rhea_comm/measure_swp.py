#!/usr/bin/env python3
'''Sweeep measurement.'''
import sys
from struct import pack
from os.path import isfile
from time    import strftime
from argparse import ArgumentParser
import numpy as np

from fpga_control  import FPGAControl
from packet_reader import read_packet_in_swp
from common import packet_size

class SwpError(Exception):
    '''Exception from sweep measurement.'''


## main
def measure_swp(fpga:FPGAControl, max_ch, f_start, f_end, f_step, fname, power, verbose=True):
    '''Do frequency sweep measurement.

    Parameters
    ----------
    fpga : FPGAControl
        FPGA controller.
    max_ch : int
        Number of DDS/DDC channels.
    f_start : float
        Start frequency in MHz.
    f_end : float
        Stop frequency in MHz.
    fname : str
        File path.
    power : int
        Number of DDSes used for each tone.
    '''

    def _vprint(*pargs, **pkwargs):
        if verbose:
            print(*pargs, **pkwargs)

    _vprint('SWEEP MEASUREMENT')
    _vprint(f'SwpPower: {power:d}/{max_ch:d}')

    fpga.init()
    fpga.iq_setting.set_read_width(1)
    cnt_finish = packet_size(1) * 10

    file_desc = open(fname, 'wb')
    fpga.tcp.clear()

    try:
        for freq in np.arange(f_start, f_end, f_step):
            _vprint(f'{freq:8.3f} MHz')
            freq_hz = int(np.floor(freq * 1e6 + 0.5))
            cnt = 0

            while True:
                fpga.dds_setting.set_freqs([freq_hz] * power + [0.] * (max_ch - power))
                fpga.dds_setting.set_amps([1]*max_ch)
                fpga.iq_setting.time_reset()
                fpga.iq_setting.iq_on()

                try:
                    while True:
                        time = read_packet_in_swp(fpga.tcp.read(packet_size(1)))
                        if time == 0:
                            break
                    break
                except KeyboardInterrupt as err:
                    print('stop timestamp reset')
                    raise KeyboardInterrupt from err
                else:
                    fpga.iq_setting.iq_off()
                    fpga.tcp.clear()
                    print('... retry read_packet_in_swp()')
                    continue


            freq_pack  = (b'\x00' * 3) if freq_hz >= 0 else (b'\xff' * 3)
            freq_pack += pack(">i", freq_hz)
            dummy_packet  = b'\xff'        # header
            dummy_packet += b'\x00' * 5    # time
            dummy_packet += freq_pack * 2  # freq
            dummy_packet += b'\xee'        # footer
            file_desc.write(dummy_packet)

            while True:
                buff = fpga.tcp.read(min(1024, cnt_finish - cnt))
                file_desc.write(buff)

                if len(buff) == 0:
                    break

                cnt += len(buff)
                if cnt >= cnt_finish:
                    break

            fpga.iq_setting.iq_off()
    except KeyboardInterrupt:
        print('stop measurement')
        fpga.iq_setting.iq_off()
    finally:
        file_desc.close()
        fpga.dac_setting.txenable_off()
        print(f'write raw data to {fname}')


def main():
    '''Parse arguments and do sweep measurement.'''
    parser = ArgumentParser()

    parser.add_argument('f_start',
                        type=float,
                        help='sweep start frequency_MHz.')

    parser.add_argument('f_end',
                        type=float,
                        help='sweep end frequency_MHz.')

    parser.add_argument('f_step',
                        type=float,
                        help='sweep start frequency_MHz.')

    parser.add_argument('-f', '--fname',
                        type=str,
                        default=None,
                        help='output filename. (default=swp_START_STOP_STEP_DATE.rawdata)')

    parser.add_argument('-p', '--power',
                        type=int,
                        default=1,
                        help='# of ch used for each comm.(<= max_ch in FPGA). default=1')

    parser.add_argument('-ip', '--ip_address',
                        type=str,
                        default='192.168.10.16',
                        help='IP-v4 address of target SiTCP. (default=192.168.10.16)')

    args = parser.parse_args()

    f_start     = args.f_start
    f_end       = args.f_end
    f_step      = args.f_step
    fname       = args.fname
    power       = args.power
    ip_address  = args.ip_address

    try:
        fpga = FPGAControl(ip_address=ip_address)
        max_ch = fpga.max_ch

        if power < 1 or power > max_ch:
            raise Exception(f'Invalid SwpPower: {power:d} / {max_ch:d}')
        if fname is None:
            fname  = 'swp'
            fname += f'_{f_start:+08.3f}MHz'
            fname += f'_{f_end:+08.3f}MHz'
            fname += f'_{f_step:+08.3f}MHz'
            fname += strftime('_%Y-%m%d-%H%M%S')
            fname += '.rawdata'
        if isfile(fname):
            raise SwpError(f'{fname!r} exists.')
    except TimeoutError:
        print('connection to FAGA failed.')
        print(ip_address, 'is invalid ip_address address.')
        sys.exit(1)
    except SwpError as err:
        print(err)
        sys.exit(1)

    measure_swp(fpga    = fpga,
                max_ch  = max_ch,
                f_start = f_start,
                f_end   = f_end,
                f_step  = f_step,
                fname   = fname,
                power   = power)



if __name__ == '__main__':
    main()
