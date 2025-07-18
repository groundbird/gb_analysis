#!/usr/bin/env python3
from time import sleep
import math

## config
cnt_step_per_sec    = 1
rate_kSPS_default   = 1 # default: 1 kSPS ~ downsample:200000
data_length_default = rate_kSPS_default * 1000 * 10 # default: 10 sec


## constant
from fpga_control import FPGAControl

## tool
def isPowerOf2(value): return value & (value-1)


## main
def output_tones(fpga:FPGAControl, MAX_CH, dds_f_MHz, power, amps=None, phases=None):
    from math import floor

    print('OUTPUT TONES')

    if amps is None:
        amps = [1.]*len(dds_f_MHz)
        for i, freq in enumerate(dds_f_MHz):
#            print(f'ch{i:03d}: {freq} MHz')
            pass
    else:
        for i, freq in enumerate(dds_f_MHz):
#            print(f'ch{i:03d}: {freq} MHz, scaled by amp')
            pass

    if phases is None:
        phases = [0.]*len(dds_f_MHz)
        for i, freq in enumerate(dds_f_MHz):
#            print(f'ch{i:03d}: {freq} MHz')
            pass
    else:
        for i, freq in enumerate(dds_f_MHz):
#            print(f'ch{i:03d}: {freq} MHz, changed by phase')
            pass


    fpga.init()
    fpga.iq_setting.set_read_width(len(dds_f_MHz))

    if power<0:
        dds_f_Hz_multi  = [freq * 1e6 for freq in dds_f_MHz]
        dds_f_Hz_multi *= floor(float(MAX_CH) / float(len(dds_f_MHz))+0.1)
        print('Power: %d*%d/%d' % (len(dds_f_MHz), MAX_CH / len(dds_f_MHz), MAX_CH))
    else:
        dds_f_Hz_multi = [freq * 1e6 for freq in dds_f_MHz] * int(power) + [0.] * (MAX_CH - int(power)*len(dds_f_MHz))
        print('Power: %d*%d/%d' % (power, len(dds_f_MHz), MAX_CH))

#    print("FREQ LIST")
#    print(dds_f_Hz_multi)
    fpga.dds_setting.set_freqs(dds_f_Hz_multi)

    if amps is not None:
        if power<0:
            amp_multi  = [amp for amp in amps]
            amp_multi *= floor(float(MAX_CH) / float(len(dds_f_MHz))+0.1)
        else:
            amp_multi = [amp for amp in amps] * int(power) + [0.] * (MAX_CH - int(power)*len(dds_f_MHz))

#        print("SCALE LIST")
#        print(amp_multi)
        fpga.dds_setting.set_amps(amp_multi)


    if phases is not None:
        if power<0:
            phase_multi = [phase for phase in phases]
            phase_multi *= floor(float(MAX_CH) / float(len(dds_f_MHz)) + 0.1)
#            print(len(phase_multi))
        else:
            phase_multi =  [phase for phase in phases] * int(power) + [0.] * (MAX_CH - int(power)*len(dds_f_MHz))

#        print("phase list")
#        print(phase_multi)
        fpga.dds_setting.set_phases(phase_multi)

## print [freq, amp, phase]
    print('INPUT list of freq, amp, phase')
    for i in range(len(dds_f_MHz)):
        print(f'ch{i:03d}: freq {dds_f_MHz[i]:7.4f}MHz, amp {amps[i]:.4f}, phase {phases[i]:7.4f}rad')

    fpga.ds_setting.set_accum(200000)

    fpga.tcp.clear()
    fpga.iq_setting.iq_on()
    _t = 0
    try:
        while True:
            buff = fpga.tcp.read(1024)
            sleep(1)
            print(f'now. {_t} s')
            _t += 1
    except Exception as e:
        print(e)
    finally:
        fpga.iq_setting.iq_off()
        fpga.dac_setting.txenable_off()
        exit(1)
    pass



if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('freqs',
                        type=float,
                        nargs='+',
                        help='list of IF target frequencys(MHz).')

    parser.add_argument('-p', '--power',
                        type=int,
                        default=1,
                        help='# of ch used for each comm.(<= max_ch in FPGA). (default: -p=1)')

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
                        dest='ip_address',
                        default='192.168.10.16',
                        help='IP_ADDRESS of target SiTCP. (default=192.168.10.16)')

    args = parser.parse_args()

    freqs = args.freqs
    amps = args.amplitude
    ip = args.ip_address

    phases = args.phase

    try:
        fpga = FPGAControl(ip_address=ip)
        MAX_CH = fpga.max_ch

        if len(freqs) > MAX_CH:
            raise Exception(f'# of channel is over max: {MAX_CH}')

        def two_div(val):
            return len(bin(val-1)) -2
        input_len = 2**two_div(len(freqs))

        if args.power > 0 and input_len > MAX_CH:
            raise Exception(f'exceeding max # of channels = {len(freqs)}*{args.power} > {MAX_CH}')

        if amps is not None and len(amps) != len(freqs):
            raise Exception(f'# of input amps should be same as # of input freqs.')
        if phases is not None and len(phases) != len(freqs):
            raise Exception(f'# of input phases must be same as # of input freqs.')

        while len(freqs) & (len(freqs)-1):
            freqs.append(0.)
            if amps is not None: amps.append(0.)
            if phases is not None: phases.append(0.)
            pass

    except TimeoutError:
        print('connection to FAGA failed.')
        print(ip, 'is invalid ip address.')
        exit(1)
        pass
    except Exception as e:
        print(e)
        exit(1)

    output_tones(fpga      = fpga,
                 MAX_CH    = MAX_CH,
                 dds_f_MHz = freqs,
                 power     = args.power,
                 amps      = amps,
                 phases     = phases)
#####
    # from sys     import argv, stderr

    # try:
    #     args        = argv[1:]
    #     arg_val     = []
    #     power       = -1
    #     amps_val    = None
    #     nch         = 0
    #     while args:
    #         if args[0] == '-p':
    #             power = int(args[1])
    #             args = args[2:]
    #         elif args[0] == '-a':
    #             print(f'AMPLITUDE SCALING IS ENABLED')
    #             if nch == 0:
    #                 raise Exception(f'ERROR:: \'-a\' option shuold be put after input frequency list')
    #             if nch+1 > len(args):
    #                 raise Exception(f'ERROR:: # of input amp should be same as # of input freqs.')
    #             amps_val = args[1:nch+1]
    #             def is_num(s):
    #                 try: float(s)
    #                 except ValueError: return False
    #                 else: return True
    #             if not all([is_num(x) for x in amps_val]):
    #                 raise Exception(f'ERROR:: non-digit input in amplitude list: {amps_val}')
    #             args = args[nch+1:]
    #         else:
    #             arg_val += [args[0]]
    #             args = args[1:]
    #             nch += 1
    #             pass
    #         pass
    #     def is_num(s):
    #         try: float(s)
    #         except ValueError: return False
    #         else: return True
    #     if not all([is_num(x) for x in arg_val]):
    #         raise Exception(f'non-digit input in frequency list: {arg_val}')
    #     if len(arg_val) > MAX_CH:
    #         raise Exception(f'too many frequencies! = {len(arg_val)} / {MAX_CH}')
    #     if power>0 and len(arg_val)*power > MAX_CH:
    #         raise Exception(f'exceeding max # of channels = {len(arg_val)}*{power} > {MAX_CH}')
    #     if len(arg_val) == 0    :
    #         raise Exception(f'no input frequencies!')
    #     def two_div(val):
    #         return len(bin(val-1)) -2
    #     input_len = 2**two_div(len(arg_val))
    #     if power >0 and power*input_len > MAX_CH:
    #         raise Exception(f'exceeding max # of channels = {input_len}*{power} > {MAX_CH}')
    #     dds_f_MHz = [float(v) for v in arg_val]
    #     if amps_val is not None and len(amps_val) != len(dds_f_MHz):
    #         raise Exception(f'# of input amp should be same as # of input freqs.')

    #     while len(dds_f_MHz) & (len(dds_f_MHz)-1):
    #         dds_f_MHz.append(0.)
    #         if amps_val is not None: amps_val.append(0.)
    # except Exception as e:
    #     print( 'Usage: measure_tod.py [(float)freq0_MHz] [(float)freq1_MHz] ...', file=stderr)
    #     print(f'               -p [(int)num] (default: {power:d})', file=stderr)
    #     print(f'                  >> amplitude (<={MAX_CH:d})', file=stderr)
    #     print(f'               -a [(float)scale0] [(float)scale1] ... (default: None)', file=stderr)
    #     print(f'                  >> amplitude scale (0.0 to 1.0). This option should be put after the frequency list.', file=stderr)
    #     print(f'                     # of amplitude scale should be same as # of input freqs.', file=stderr)
    #     print(e)
    #     exit(1)

    # output_tones(dds_f_MHz, power, amps=amps_val)
