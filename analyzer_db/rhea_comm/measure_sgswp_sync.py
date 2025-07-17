#!/usr/bin/env python3
import re

from fpga_control import fpga_control

class FUnit(object):
    mHz, Hz, kHz, MHz, GHz = 1e-3, 1.0, 1e+3, 1e+6, 1e+9

def freq_fmt_Hz(f_num, hz_scale):
    return round(f_num*hz_scale)

def fstr_to_int(f_str, unit_ret=None):
    result = re.search(r'\d+(\.\d+)?', f_str)    
    f_num_str = result.group()
    unit_str = f_str[result.end():]
    if unit_ret is None: unit_ret = 'Hz'
    unit = getattr(FUnit, unit_str) / getattr(FUnit, unit_ret)
    f_ret = float(f_num_str) * unit
    return round(f_ret)

def measure_sgswp(fpga,
                  MAX_CH,
                  dds_f_MHz,
                  data_length,
                  rate_kSPS,
                  power,
                  fname,
                  f_start_MHz,
                  f_stop_MHz,
                  f_step_kHz,
                  amps=None):
    from math import floor
    from struct  import pack

    print('SG_SWP MEASUREMENT')
    if amps is None:
        amps = [1.]*len(dds_f_MHz)
        for i, freq in enumerate(dds_f_MHz):
            print(f'ch{i:03d}: {freq} MHz')
            pass
    else:
        for i, freq in enumerate(dds_f_MHz):
            print(f'ch{i:03d}: {freq} MHz, scaled by amp')
            pass

    fpga.init()
    fpga.iq.set_read_width(len(dds_f_MHz))

    if power<0:
        dds_f_Hz_multi  = [freq * 1e6 for freq in dds_f_MHz]
        dds_f_Hz_multi *= floor(float(MAX_CH) / float(len(dds_f_MHz))+0.1)
        print('Power: %d*%d/%d' % (len(dds_f_MHz), MAX_CH / len(dds_f_MHz), MAX_CH))
    else:
        dds_f_Hz_multi = [freq * 1e6 for freq in dds_f_MHz] * int(power) + [0.] * (MAX_CH - int(power)*len(dds_f_MHz))
        print('Power: %d*%d/%d' % (power, len(dds_f_MHz), MAX_CH))

    print("INPUT FREQ LIST")
    print(dds_f_Hz_multi)
    fpga.dds.set_freqs(dds_f_Hz_multi)

    if amps is not None:
        if power<0:
            amp_multi  = [amp for amp in amps]
            amp_multi *= floor(float(MAX_CH) / float(len(dds_f_MHz))+0.1)
        else:
            amp_multi = [amp for amp in amps] * int(power) + [0.] * (MAX_CH - int(power)*len(dds_f_MHz))

        print("INPUT SCALE LIST")
        print(amp_multi)
        fpga.dds.set_amps(amp_multi)

    fpga.ds.set_rate(floor(200000 / rate_kSPS+0.5))

    f = open(fname, 'wb')

    ## write header
    dummy_packet  = b'\xaa' # header
    dummy_packet += b'\x00' + pack('>I', rate_kSPS * 1000) # timestamp
    for freq in dds_f_MHz:
        freq_Hz = int(floor(freq * 1e6 + 0.5))
        freq_packet  = (b'\x00' * 3) if freq_Hz >= 0 else (b'\xff' * 3)
        freq_packet += pack('>i', freq_Hz)
        swp_packet  = b'\x00'
        swp_packet += pack('>H', f_start_MHz)
        swp_packet += pack('>H', f_stop_MHz)
        swp_packet += pack('>H', f_step_kHz)
        dummy_packet += freq_packet
        dummy_packet += swp_packet
        pass
    dummy_packet += b'\xee' # footer
    f.write(dummy_packet)

    cnt = 0
    packet_size = 7 + 7 * 2 * len(dds_f_MHz)
    cnt_finish = packet_size * data_length
    cnt_step = packet_size * rate_kSPS * 1000
    cnt_print = cnt_step

    fpga.tcp.clear()
    fpga.iq.iq_on()

    try:
        while True:
            buff = fpga.tcp.read(min(1024, cnt_finish - cnt))
            f.write(buff)
            if not len(buff): break
            cnt += len(buff)
            while cnt >= cnt_print:
                print(cnt_print / packet_size)
                cnt_print += cnt_step
                pass
            if cnt >= cnt_finish: break
            pass
    except KeyboardInterrupt:
        print('stop measurement')
    finally:
        f.close()
        fpga.iq.iq_off()
        fpga.dac.txenable_off()
        print(f'write raw data to {fname}')
        pass


if __name__ == '__main__':
    from os.path  import isfile
    from argparse import ArgumentParser
    from time import strftime

    parser = ArgumentParser()

    parser.add_argument('freqs',
                        type=float,
                        nargs='+',
                        help='list of target frequencys(MHz).')
    
    parser.add_argument('-f', '--fname',
                        type=str,
                        default=None,
                        help='user setting filename')
    
    parser.add_argument('-r', '--rate',
                        type=int,
                        default=10,
                        help=f'rate(kSPS). default: {10}')
    
    parser.add_argument('-p', '--power',
                        type=int,
                        default=-1,
                        help='number')
    
    parser.add_argument('-a', '--amps',
                        type=float,
                        nargs='+',
                        default=None,
                        help='''list of float, amplitude scale(0.0 to 1.0).
                        # of amplitude scale should be same as # of input freqs.''')

    parser.add_argument('-fc', '--f_center',
                        type=str,
                        dest='f_center',
                        default='4GHz',
                        help='center frequency of sweep with unit.')

    parser.add_argument('-fw', '--f_width',
                        type=str,
                        dest='f_width',
                        default='2MHz',
                        help='sweeping frequency half-width with unit.')

    parser.add_argument('-fs', '--f_step',
                        type=str,
                        dest='f_step',
                        default='10kHz',
                        help='sweep step frequency with unit.')

    parser.add_argument('-d', '--dwell',
                        type=int,
                        dest='dwell',
                        default=1000,
                        help='sweep dwell time in (us).')

    parser.add_argument('-nr', '--n_run',
                        type=int,
                        dest='n_run',
                        default=1,
                        help='# of sweep run.')

    parser.add_argument('--mon',
                        type=bool,
                        dest='mon',
                        default=False,
                        help='if use as sweep monitor, --mon True.')

    parser.add_argument('--ip',
                        type=str,
                        dest='ip',
                        default='192.168.10.16',
                        help='IP_ADDRESS of target SiTCP')

    # parser.add_argument('--fastswp',
    #                     action='store_true',
    #                     help='flag of using fast sweep mode.')

    args = parser.parse_args()
    freqs = args.freqs
    amps = args.amps
    fname = args.fname

    fpga = fpga_control(ip_address=args.ip)
    MAX_CH = fpga.MAX_CH

    try:
        if len(freqs) > MAX_CH :
            raise Exception(f'# of channel is over max: {MAX_CH}')

        def two_div(val):
            return len(bin(val-1)) -2
        input_len = 2**two_div(len(freqs))

        if args.power > 0 and input_len > MAX_CH:
            raise Exception(f'exceeding max # of channels = {len(freqs)}*{args.power} > {MAX_CH}')

        if 200000 % args.rate != 0:
            raise Exception(f'sampling rate [kSPS] should be a divisor of 2000000: input rate = {args.rate} kSPS')

        if fname is None:
            # if args.fastswp is False:
            #     fname  = 'sgswp'
            # else:
            #     fname  = 'sgfastswp'
            fname  = 'sgswp'
            for f in args.freqs:
                fname += f'_{f:+08.3f}MHz'
                pass
            fname += f'_{args.rate:04d}kSPS'
            fname += f'_{args.f_center}center_{args.f_width}width_{args.f_step}_step'
            fname += strftime('_%Y-%m%d-%H%M%S')
            fname += '.rawdata'
        if args.mon is True:
            pass
        elif isfile(fname):
            raise Exception(f'{fname} is existed.')

        if amps is not None and len(amps) != len(freqs):
            raise Exception(f'# of input amps should be same as # of input freqs.') 

        while len(freqs) & (len(freqs)-1):
            freqs.append(0.)
            if amps is not None: amps.append(0.)

    except Exception as e:
        print(e)
        exit(1)

    # if args.fastswp is False:
    #     fineness = fstr_to_int(args.f_step, 'kHz')
    # else:
    #     fineness = fstr_to_int(args.f_step + 'kHz', 'kHz')
    #     pass

    f_center_MHz = fstr_to_int(args.f_center, 'MHz')
    f_width_MHz = fstr_to_int(args.f_width, 'MHz')
    f_start_MHz = f_center_MHz - f_width_MHz
    f_stop_MHz  = f_center_MHz + f_width_MHz
    f_step_kHz = fstr_to_int(args.f_step, 'kHz')
    length = (f_width_MHz * 2 / f_step_kHz * (args.dwell+200) * args.n_run * 1e-3 + 2) * args.rate * 1000

    measure_sgswp(fpga        = fpga,
                  MAX_CH      = MAX_CH,
                  dds_f_MHz   = freqs,
                  data_length = length,
                  rate_kSPS   = args.rate,
                  power       = args.power,
                  fname       = fname,
                  f_start_MHz = f_start_MHz,
                  f_stop_MHz  = f_stop_MHz,
                  f_step_kHz  = f_step_kHz,
                  amps        = amps)
