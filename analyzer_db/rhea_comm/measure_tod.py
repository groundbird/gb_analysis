#!/usr/bin/env python3
'''Measure time-ordered data (TOD).'''
from math import floor
from struct  import pack
from time    import strftime
from os.path import isfile
from argparse import ArgumentParser
import sys

from fpga_control import FPGAControl
from common import two_div, packet_size

## config
CNT_STEP_PER_SEC    = 1
RATE_KSPS_DEFAULT   = 1
DATA_LENGTH_DEFAULT = RATE_KSPS_DEFAULT * 1000 * 10 # default: 10 sec
READ_RATE = 0.1 # sec

class TODError(Exception):
    '''Raised when error happens during TOD measurement.'''

## main
def measure_tod(fpga:FPGAControl, max_ch, dds_f_megahz, data_length,
                rate_ksps, power, fname, amps=None, phases=None, verbose=True, swap_dac=True, swap_adc=True):
    '''Measure time-ordered data.

    Parameters
    ----------
    fpga : FPGAControl
        FPGA controller.
    max_ch : int
        number of DDS channels in FPGA.
    dds_f_megahz : list of float
        DDS frequency in MHz.
    data_length : int
        Data acquisition length.
    rate_ksps : int
        Data sampling rate in kHz.
    power : int
        Number of DDSes used for each tone.
    fname : str
        File name.
    amps : list of float
        Amplitude of tones.
    phases : list of float
        DDS initial phases in radian.
    swap_dac : boolean, optional
        Whether I and Q for DAC are swapped or not.
    swap_adc : boolean, optional
        Whether I and Q for ADC are swapped or not.
    '''
    def _vprint(*pargs, **pkwargs):
        if verbose:
            print(*pargs, **pkwargs)

    _vprint('TOD MESUREMENT.')

    if amps is None:
        amps = [1.]*len(dds_f_megahz)

    if phases is None:
        phases = [0.]*len(dds_f_megahz)
    fpga.init(swap_dac=swap_dac, swap_adc=swap_adc)
    fpga.iq_setting.set_read_width(len(dds_f_megahz))

    if power < 0:
        dds_f_hz_multi  = [freq * 1e6 for freq in dds_f_megahz]
        dds_f_hz_multi *= floor(float(max_ch) / float(len(dds_f_megahz))+0.1)
    else:
        dds_f_hz_multi = [freq * 1e6 for freq in dds_f_megahz] * int(power)
        dds_f_hz_multi += [0.] * (max_ch - len(dds_f_hz_multi))

    fpga.dds_setting.set_freqs(dds_f_hz_multi)

    if amps is not None:
        if power < 0:
            amp_multi  = list(amps)
            amp_multi *= floor(float(max_ch) / float(len(dds_f_megahz))+0.1)
        else:
            amp_multi = list(amps) * int(power)
            amp_multi += [0.] * (max_ch - len(amp_multi))

        fpga.dds_setting.set_amps(amp_multi)

    if phases is not None:
        if power<0:
            phase_multi = list(phases)
            phase_multi *= floor(float(max_ch) / float(len(dds_f_megahz)) + 0.1)
        else:
            phase_multi = list(phases) * int(power)
            phase_multi += [0.] * (max_ch - len(phase_multi))

        fpga.dds_setting.set_phases(phase_multi)


    _vprint('INPUT list of freq, amp, phase')
    for i, freq in enumerate(dds_f_megahz):
        _vprint(f'ch{i:03d}: freq {freq:7.4f}MHz, amp {amps[i]:.4f}, phase {phases[i]:7.4f}rad ')

    fpga.ds_setting.set_accum(floor(200000 / rate_ksps+0.5))

    file_desc = open(fname, 'wb')

    ## write header
    dummy_packet  = b'\xff'
    dummy_packet += b'\x00' + pack('>I', rate_ksps * 1000)

    for freq in dds_f_megahz:
        freq_hz = int(floor(freq * 1e6 + 0.5))
        freq_packet = (b'\x00' * 3) if freq_hz >= 0 else (b'\xff' * 3)
        freq_packet += pack('>i', freq_hz)
        dummy_packet += freq_packet * 2

    dummy_packet += b'\xee'
    file_desc.write(dummy_packet)

    cnt = 0
    psize = packet_size(len(dds_f_megahz))
    cnt_finish = psize * data_length
    cnt_step = psize * rate_ksps * 1000 * CNT_STEP_PER_SEC
    cnt_print = cnt_step

    # Buffer size calculation.
    buffsize = psize * rate_ksps * 1000 * READ_RATE
    buffsize = 1024 if buffsize < 1024 else int(buffsize)

    fpga.tcp.clear()
    fpga.iq_setting.iq_on()

    try:
        while True:
            buff = fpga.tcp.read(min(buffsize, cnt_finish - cnt))
            file_desc.write(buff)

            if len(buff) == 0:
                break

            cnt += len(buff)

            while cnt >= cnt_print:
                _vprint(cnt_print / psize)
                cnt_print += cnt_step

            if cnt >= cnt_finish:
                break

    except KeyboardInterrupt:
        print('stop measurement')
    finally:
        file_desc.close()
        fpga.iq_setting.iq_off()
        fpga.dac_setting.txenable_off()
        _vprint(f'write raw data to {fname}')


def main():
    '''TOD measurement'''
    parser = ArgumentParser()

    parser.add_argument('freqs',
                        type=float,
                        nargs='+',
                        help='list of target frequency_MHz.')

    parser.add_argument('-f', '--fname',
                        type=str,
                        default=None,
                        help='output filename. (default=tod_INPUTFREQ_RATE_DATE.rawdata)')

    parser.add_argument('-l', '--length',
                        type=int,
                        default=DATA_LENGTH_DEFAULT,
                        help=f'data length. (default={DATA_LENGTH_DEFAULT})')

    parser.add_argument('-r', '--rate',
                        type=int,
                        default=RATE_KSPS_DEFAULT,
                        help=f'Sampling rate in kSPS. (default={RATE_KSPS_DEFAULT}kSPS)')

    parser.add_argument('-p', '--power',
                        type=int,
                        default=1,
                        help='# of ch used for each comm.(<= max_ch in FPGA). default=1')

    parser.add_argument('--amplitude',
                        type=float,
                        nargs='+',
                        default=None,
                        help='''list of amplitude scale (0 to 1.0).
                        # of amplitude scale must be same as # of input freqs.''')

    parser.add_argument('--phase',
                        type=float,
                        nargs='+',
                        default=None,
                        help='''list of phase scale[rad].
                        # of phase scale must be same as # of input freqs''')

    parser.add_argument('-ip', '--ip_address',
                        type=str,
                        default='192.168.10.16',
                        help='IP-v4 address of target SiTCP. (default=192.168.10.16)')

    args = parser.parse_args()

    freqs       = args.freqs
    fname       = args.fname
    data_length = args.length
    rate_ksps   = args.rate
    power       = args.power
    amps        = args.amplitude
    phases      = args.phase
    ip_address  = args.ip_address

    try:
        fpga = FPGAControl(ip_address=ip_address)
        max_ch = fpga.max_ch

        if len(freqs) > max_ch:
            raise TODError(f'too many frequencies! = {len(freqs)} / {max_ch}')
        if power>0 and len(freqs)*power > max_ch:
            raise TODError(f'exceeding max # of channels = {len(freqs)}*{power} > {max_ch}')

        input_len = 2**two_div(len(freqs))
        if power >0 and power*input_len > max_ch:
            raise TODError(f'exceeding max # of channels = {input_len}*{power} > {max_ch}')

        if 200000 % rate_ksps != 0:
            raise TODError('sampling rate [kSPS] should be a divisor of 2000000:'
                            + f'input rate = {rate_ksps} kSPS')

        dds_f_megahz = [float(v) for v in freqs]

        if fname is None:
            fname  = 'tod'
            for freq in dds_f_megahz:
                fname += f'_{freq:+08.3f}MHz'

            fname += f'_{rate_ksps:04d}kSPS'
            fname += strftime('_%Y-%m%d-%H%M%S')
            fname += '.rawdata'

        if isfile(fname):
            raise Exception(f'{fname} exists.')

        if amps is not None and len(amps) != len(freqs):
            raise Exception('# of input amp must be same as # of input freqs.')

        if phases is not None and len(phases) != len(freqs):
            raise Exception('# of input phase must be same as # of input freqs.')

    except TimeoutError:
        print('connection to FPGA failed.')
        print(ip_address, 'is invalid ip address.')
        sys.exit(1)
    except TODError as err:
        print(err)
        sys.exit(1)


    while len(dds_f_megahz) & (len(dds_f_megahz)-1):
        dds_f_megahz.append(0.)
        if amps is not None:
            amps.append(0.)
        if phases is not None:
            phases.append(0.)

    #fpga.init()
    measure_tod(fpga        = fpga,
                max_ch      = max_ch,
                dds_f_megahz   = dds_f_megahz,
                data_length = data_length,
                rate_ksps   = rate_ksps,
                power       = power,
                fname       = fname,
                amps        = amps,
                phases       = phases)


if __name__ == '__main__':
    main()
