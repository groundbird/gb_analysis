#!/usr/bin/env python3
import re
from time import sleep
from sg_sweep_sync import Sweep

## constant
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
                  amps=None,
                  swp=None):
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
    sleep(0.001)
    tmp = swp.send_command()

    try:
        while True:
            buff = fpga.tcp.read(min(1024, cnt_finish - cnt))
            f.write(buff)
            if not len(buff): break
            cnt += len(buff)
            while cnt >= cnt_print:
                print(round(cnt_print/cnt_step), '/', round(cnt_finish/cnt_step))
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
    from time import strftime
    from argparse import ArgumentParser
    import serial.tools.list_ports as ports

    parser = ArgumentParser()

    parser.add_argument('freqs',
                        type=float,
                        nargs='+',
                        help='list of IF target frequencys(MHz).')

    parser.add_argument('-f', '--fname',
                        type=str,
                        default=None,
                        help='output filename. (default=sgswp_IFFREQ_RATE_CENTER_WIDTH_STEP_DATE.rawdata)')

    parser.add_argument('-r', '--rate',
                        type=int,
                        default=10,
                        help='sample rate of measurement (kSPS). (default=10kSPS)')

    parser.add_argument('-p', '--power',
                        type=int,
                        default=-1,
                        help='# of ch used for each comm.(<= max_ch in FPGA). (default=MAX_CH/n_comm)')

    parser.add_argument('-a', '--amplitude',
                        type=float,
                        nargs='+',
                        default=None,
                        help='''list of amplitude scale (0 to 1.0).
                        # of amplitude scale must be same as # of input freqs.''')

    parser.add_argument('-fc', '--f_center',
                        type=str,
                        dest='f_center',
                        default='4GHz',
                        help='center frequency of sweep with unit. (default=4GHz)')

    parser.add_argument('-fw', '--f_width',
                        type=str,
                        dest='f_width',
                        default='2MHz',
                        help='sweeping frequency half-width with unit (sweep range = fcenter +/- fwidth). (default=2MHz)')

    parser.add_argument('-fs', '--f_step',
                        type=str,
                        dest='f_step',
                        default='10kHz',
                        help='sweep step frequency with unit. (default=10kHz)')

    parser.add_argument('-d', '--dwell',
                        type=int,
                        dest='dwell',
                        default=1000,
                        help='sweep dwell time in (us). (default=1000us=1ms)')

    parser.add_argument('-nr', '--n_run',
                        type=int,
                        dest='n_run',
                        default=1,
                        help='# of sweep run. (default=1)')

    parser.add_argument('--port',
                        type=str,
                        dest='port',
                        default='/dev/ttyS5',
                        help='''path of COM port for SG-control.
                        for WSL: /dev/ttyS_, for mac: /dev/tts.usbmodem_, for linux: /dev/ttyACM_.
                        (default=/dev/ttyS5)''')

    parser.add_argument('-ip', '--ip_address',
                        type=str,
                        dest='ip_address',
                        default='192.168.10.16',
                        help='IP_ADDRESS of target SiTCP. (default=192.168.10.16)')

    parser.add_argument('--mon',
                        type=bool,
                        dest='mon',
                        default=False,
                        help='if use as sweep monitor, --mon True. (default=False)')

    args = parser.parse_args()
    freqs = args.freqs
    amps = args.amplitude
    fname = args.fname
    ip = args.ip_address

    try:
        fpga = fpga_control(ip_address=ip)
        MAX_CH = fpga.MAX_CH

        if len(freqs) > MAX_CH:
            raise Exception(f'# of channel is over max: {MAX_CH}')

        def two_div(val):
            return len(bin(val-1)) -2
        input_len = 2**two_div(len(freqs))

        if args.power > 0 and input_len > MAX_CH:
            raise Exception(f'exceeding max # of channels = {len(freqs)}*{args.power} > {MAX_CH}')

        if 200000 % args.rate != 0:
            raise Exception(f'sampling rate [kSPS] should be a divisor of 2000000: input rate = {args.rate} kSPS')

        if fname is None:
            fname  = 'sgswp'
            for f in args.freqs:
                fname += f'_{f:+08.3f}MHz'
                pass
            fname += f'_{args.rate:04d}kSPS'
            fname += f'_{args.f_center}center_{args.f_width}width_{args.f_step}_step'
            fname += strftime('_%Y-%m%d-%H%M%S')
            fname += '.rawdata'
        if (args.mon is False) and isfile(fname):
            raise Exception(f'{fname} is existed.')

        if amps is not None and len(amps) != len(freqs):
            raise Exception(f'# of input amps should be same as # of input freqs.')

        while len(freqs) & (len(freqs)-1):
            freqs.append(0.)
            if amps is not None: amps.append(0.)

        coms = ports.comports()
        com_str = []
        for com in coms:
            com_str.append(com.device)
            pass
        if args.port not in com_str:
            raise Exception(args.port, 'is invalid port path.')

    except TimeoutError:
        print('connection to FAGA failed.')
        print(ip, 'is invalid ip address.')
        exit(1)
        pass
    except Exception as e:
        print(e)
        exit(1)

    f_center_MHz = fstr_to_int(args.f_center, 'MHz')
    f_width_MHz = fstr_to_int(args.f_width, 'MHz')
    f_start_MHz = f_center_MHz - f_width_MHz
    f_stop_MHz  = f_center_MHz + f_width_MHz
    f_step_kHz = fstr_to_int(args.f_step, 'kHz')
    t_step_ms = args.dwell * 1e-3 + 0.12 # locked time ~ 110us
    length = (f_width_MHz * 2 / f_step_kHz * args.n_run * t_step_ms + 0.001) * 1000 * args.rate

    sweep = Sweep(port   = args.port,
                  mode   = 'normal',
                  start  = f_start_MHz*1e6,
                  stop   = f_stop_MHz*1e6,
                  step   = f_step_kHz*1e3,
                  points = 0,
                  dwell  = args.dwell,
                  run    = args.n_run+1)

    print('Target SG Serial # :', sweep.get_id().serial_number)

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
                  amps        = amps,
                  swp         = sweep)

    print(f'SG: f_center={f_center_MHz}MHz, f_width={f_width_MHz}MHz, f_step={f_step_kHz}kHz')
    sweep.set_frequency(args.f_center)
