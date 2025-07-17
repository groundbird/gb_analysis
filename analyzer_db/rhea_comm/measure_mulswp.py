#!/usr/bin/env python3
'''Perform multi-channel sweep.'''
from os.path import isfile
from time    import strftime
from argparse import ArgumentParser
import sys
import numpy as np

from fpga_control  import FPGAControl
from packet_reader import read_packet_in_swp
from common import two_div, packet_size

class MulswpError(Exception):
    '''Exception raised by mulsweep function.'''

## main
def measure_mulswp(fpga:FPGAControl, max_ch, dds_f_megahz, width, step, fname,
                   power, amps=None, phases=None, f_off=None, verbose=True, swap_dac=True, swap_adc=True):
    '''Perform multi-channel sweep.

    Parameters
    ----------
    fpga : FPGAControl
        FPGA controller.
    max_ch : int
        Number of DDS channels in the firmware.
    dds_f_megahz : list of float
        Tone frequency list in MHz.
    width : float
        Sweep width in MHz.
    step : float
        Step frequency in MHz.
    fname : str
        File path.
    power : int
        Number of DDSes to be used for a tone.
    amps : list of float, optional
        Amplitude list (maximum: 1).
    phases : list of float, optional
        Initial phase list in radians.
    f_off : float, optional
        Frequency of off-resonance measurement in MHz.
    swap_dac : boolean, optional
        Whether I and Q for DAC are swapped or not.
    swap_adc : boolean, optional
        Whether I and Q for ADC are swapped or not.
    '''
    def _vprint(*pargs, **pkwargs):
        if verbose:
            print(*pargs, **pkwargs)

    n_tones = len(dds_f_megahz)

    if power < 0:
        num_list = np.floor(float(max_ch) / float(len(dds_f_megahz))+0.1)
    else:
        num_list = int(power)

    # Amplitude processing
    if amps is None:
        amps = [1.]*len(dds_f_megahz)

    if power < 0:
        amp_multi  = list(amps)
        amp_multi *= num_list
    else:
        amp_multi = list(amps) * num_list
        amp_multi += [0.] * (max_ch - len(amp_multi))

    # Phase processing
    if phases is None:
        phases = [0.]*len(dds_f_megahz)

    if power < 0:
        phase_multi = list(phases)
        phase_multi *= num_list
    else:
        phase_multi =  list(phases) * num_list
        phase_multi += [0.] * (max_ch - len(phase_multi))

    _vprint('INPUT list of freq, amp, phase')
    for i, (freq, amp, phase) in enumerate(zip(dds_f_megahz, amp_multi, phase_multi)):
        _vprint(f'ch{i:03d}: freq {dds_f_megahz[i]:7.4f}MHz, amp {amp:.4f}, phase {phase:7.4f}rad')

    input_len = 2**two_div(n_tones)
    if f_off is not None:
        input_len *= 2

    psize = packet_size(input_len)
    cnt_finish = 10 * psize

    _vprint('MULTI-SWEEP MEASUREMENT')
    _vprint(f'SwpPower: {power:d}*{input_len}/{max_ch:d}')

    #fpga.init()
    fpga.init(swap_dac=swap_dac, swap_adc=swap_adc)
    print(f'Input len: {input_len}')
    fpga.iq_setting.set_read_width(input_len)

    file_desc = open(fname, 'wb')
    fpga.tcp.clear()

    dfreqs = np.arange(-width/2, width/2, step)

    try:
        for dfreq in dfreqs:
            _vprint(f'{dfreq:8.3f} MHz')
            freq_hzs = [int(np.floor((freq+dfreq) * 1e6 + 0.5)) for freq in dds_f_megahz]

            if power < 0:
                input_freqs = freq_hzs * num_list
            else:
                input_freqs = freq_hzs * num_list
                input_freqs += [0]*(max_ch - len(input_freqs))

            cnt = 0

            while True:
                fpga.dds_setting.set_freqs(input_freqs)
                fpga.dds_setting.set_amps(amp_multi)
                fpga.dds_setting.set_phases(phase_multi)
                fpga.iq_setting.time_reset()
                fpga.iq_setting.iq_on()

                try:
                    while True:
                        time = read_packet_in_swp(fpga.tcp.read(psize))
                        if time == 0:
                            break
                    break
                except KeyboardInterrupt as err:
                    _vprint('stop reset timestamp')
                    raise KeyboardInterrupt from err
                else:
                    fpga.iq_setting.iq_off()
                    fpga.tcp.clear()
                    _vprint('... retry read_packet_in_swp()')

            freq_packs = [int(freq_hz).to_bytes(7, byteorder='big', signed=True) \
                          for freq_hz in freq_hzs]

            dummy_packet  = b'\xff'        # header
            dummy_packet += b'\x00' * 5    # time

            for freq_pack in freq_packs:
                dummy_packet += freq_pack * 2     # freq

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
        _vprint('stop measurement')
        fpga.iq_setting.iq_off()
    finally:
        file_desc.close()
        fpga.dac_setting.txenable_off()
        print(f'write raw data to {fname}')


def main():
    '''Parse arguments and perform multi-tone sweep measurement.'''
    parser = ArgumentParser()

    parser.add_argument('fcenters',
                        type=float,
                        nargs='+',
                        help='list of sweep center frequency_MHz.')

    parser.add_argument('-f', '--fname',
                        type=str,
                        default=None,
                        help='output filename. (default=mulswp_WIDTH_STEP_INPUTFREQ_DATE.rawdata)')

    parser.add_argument('-p', '--power',
                        type=int,
                        default=1,
                        help='# of ch used for each comm.(<= max_ch in FPGA). default=1')

    parser.add_argument('-w', '--width',
                        type=float,
                        default=3.0,
                        help='frequency width in MHz'+
                        '(frequency range = freq +/- width/2). default=3.0')


    parser.add_argument('-s', '--step',
                        type=float,
                        default=0.01,
                        help='frequency step in MHz. default=0.01')

    parser.add_argument('-o', '--offreso_tone',
                        type=float,
                        default=None,
                        help='off resonance frequency fixed to same value during sweeping.'
                              + ' default=None')

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

    dds_f_megahz = args.fcenters
    fname      = args.fname
    power      = args.power
    width      = args.width
    step       = args.step
    f_off      = args.offreso_tone
    amps       = args.amplitude
    phases     = args.phase
    ip_address = args.ip_address

    try:
        fpga = FPGAControl(ip_address=ip_address)
        max_ch = fpga.max_ch

        input_len = 2**two_div(len(dds_f_megahz)+(1 if f_off is not None else 0))
        if power < 1 or power*input_len > max_ch:
            raise MulswpError(f'exceeding max # of channels = {input_len}*{power} > {max_ch}')
        if fname is None:
            fname  = 'mulswp'
            fname += f'_{width:+08.3f}MHzWidth'
            fname += f'_{step:+08.3f}MHzStep'
            for f_cen in dds_f_megahz:
                fname += f'_{f_cen:+08.3f}MHz'
            if f_off is not None:
                fname += f'_{f_off:+08.3f}MHzOffTone'
            fname += strftime('_%Y-%m%d-%H%M%S')
            fname += '.rawdata'
        if isfile(fname):
            raise MulswpError(f'{fname} exists.')
        if amps is not None and len(amps) != len(dds_f_megahz):
            raise MulswpError('# of input amp must be same as # of input freqs.')
        if phases is not None and len(phases) != len(dds_f_megahz):
            raise MulswpError('# of input phase must be same as # of input freqs.')

    except TimeoutError:
        print('connection to FAGA failed.')
        print(ip_address, 'is invalid ip address.')
        sys.exit(1)
    except MulswpError as err:
        print(err)
        sys.exit(1)

    while len(dds_f_megahz) & (len(dds_f_megahz)-1):
        dds_f_megahz.append(0.)
        if amps is not None:
            amps.append(0.)
        if phases is not None:
            phases.append(0.)

    measure_mulswp(fpga      = fpga,
                   max_ch    = max_ch,
                   dds_f_megahz = dds_f_megahz,
                   width     = width,
                   step      = step,
                   fname     = fname,
                   power     = power,
                   amps      = amps,
                   phases    = phases,
                   f_off     = f_off)


if __name__ == '__main__':
    main()
