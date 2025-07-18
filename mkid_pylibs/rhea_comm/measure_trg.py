#!/usr/bin/env python3
'''Trigger measurement.'''
from argparse import ArgumentParser
from time import strftime, sleep
from pathlib import Path

from math    import floor
from struct  import pack
from datetime import datetime, timezone, timedelta
from numpy   import mean, std

from fpga_control import FPGAControl
from rhea_pkg import RATE_KSPS_DEFAULT, IP_ADDRESS_DEFAULT
from packet_reader import read_iq_packet
from tone_conf import ToneConf
from common import packet_size

PREMEAN_LEN_DEFAULT = 10000
DATA_LEN_DEFAULT = 1024
THR_SIGMA_DEFAULT = 3.
THR_COUNT_DEFAULT = 1
TRIG_POS_DEFAULT = 100

JST = timezone(timedelta(hours=+9), 'JST')

class TrgError(Exception):
    '''Exception raised in trigger measurement.'''


## main
def measure_trg(fpga:FPGAControl, tone_conf:ToneConf, data_length,
                thre_sigma, thre_count, rate_ksps,
                trig_pos, pre_length, fname, verbose=True,
                end=None):
    '''Perform trigger measurement.
    '''

    def _vprint(*pargs, **pkwargs):
        if verbose:
            print(*pargs, **pkwargs)

    _vprint('TRIGGER MEASUREMENT')

    dds_f_megahz = [freq/1e6 for freq in tone_conf.freq_if]

    for i, freq in enumerate(dds_f_megahz):
        _vprint(f'ch{i:03d}: {freq} MHz')

    fpga.init()
    fpga.iq_setting.set_read_width(tone_conf.n_tone)

    fpga.dds_setting.configure(tone_conf)

    fpga.ds_setting.set_accum(floor(200000 / rate_ksps+0.1))
    fpga.trg_setting.set_trig_pos(trig_pos)

    ## make dummpy packet
    dummy_packet  = b'\xff'
    dummy_packet += b'\x00' + pack('>I', rate_ksps * 1000)
    dummy_packet += tone_conf.freq_repr
    dummy_packet += b'\xee'

    ## pre-measurement
    cnt = 0
    psize = packet_size(tone_conf.n_tone)
    cnt_finish = psize * pre_length

    _vprint('pre-measurement')

    data_stock = [0 for i in range(pre_length)]
    cnt = 0

    readlen = psize*pre_length
    
    errcnt = 0
    while True:
        fpga.tcp.clear()
        fpga.iq_setting.iq_on()
        if errcnt > 10:
            raise Exception('too many error')
        try:
            rawdata = fpga.tcp.read(readlen)
            if len(rawdata) == readlen:
                break
        except Exception as err:
            print(err)
            errcnt += 1

    for i in range(pre_length):
        cnt = i
        _, data, _, _ = read_iq_packet(rawdata[psize*i:psize*(i+1)])
        data_stock[i] = data

    fpga.iq_setting.iq_off()

    _vprint(f'read {cnt:d} events')
    data_stock = list(zip(*data_stock[0:cnt]))
    mean_set  = [mean(d) for d in data_stock]
    sigma_set = [std(d)  for d in data_stock]
    _vprint('min:', end=' ')
    for _m, _s in zip(mean_set, sigma_set):
        _vprint(int(_m - thre_sigma * _s), end=' ')

    _vprint()
    _vprint('max:', end=' ')
    for _m, _s in zip(mean_set, sigma_set):
        _vprint(int(_m + thre_sigma * _s), end=' ')
    _vprint()

    ## set threshold
    for i in range(floor(len(mean_set)/2)):
        i_mean,  q_mean  = mean_set [2*i:2*i+2]
        i_sigma, q_sigma = sigma_set[2*i:2*i+2]
        i_min = int(i_mean - thre_sigma * i_sigma)
        i_max = int(i_mean + thre_sigma * i_sigma)
        q_min = int(q_mean - thre_sigma * q_sigma)
        q_max = int(q_mean + thre_sigma * q_sigma)
        fpga.trg_setting.set_threshold(i, [i_min, i_max], [q_min, q_max], trg_reset=False)

    fpga.trg_setting.set_thre_count(thre_count, trg_reset=True)

    ## get trigger event
    _vprint('main-measurement')
    fpga.tcp.clear()
    fpga.trg_setting.start()

    _vprint('wait trigger')
    cnt = 0
    while fpga.trg_setting.state() == 1:
        cnt += 1
        if cnt % 10 == 0:
            _vprint(cnt)
        else:
            _vprint('.', end='')

        if end is not None:
            if end < datetime.now(tz=JST):
                fpga.iq_setting.iq_off()
                fpga.dac_setting.txenable_off()
                fpga.tcp.clear()
                return

        sleep(1)

    _vprint('triggerd')

    file_desc = open(fname, 'wb')
    file_desc.write(dummy_packet)

    cnt = 0
    cnt_finish = psize * data_length

    while True:
        buff = fpga.tcp.read(min(1024, cnt_finish - cnt))
        file_desc.write(buff)
        if len(buff) == 0:
            break
        cnt += len(buff)
        if cnt >= cnt_finish:
            break

    file_desc.close()
    fpga.iq_setting.iq_off()
    fpga.dac_setting.txenable_off()

    _vprint(f'Write raw data to {fname}')


def main():
    '''Parse arguments and do trigger measurement.'''
    parser = ArgumentParser()

    parser.add_argument('freqs',
                        type=float,
                        nargs='+',
                        help='list of tone frequencies in MHz.')

    parser.add_argument('-f', '--fname',
                        type=str,
                        default=None,
                        help='output filename. (default=mulswp_WIDTH_STEP_INPUTFREQ_DATE.rawdata)')

    parser.add_argument('-t', '--threshold',
                        type=float,
                        default=THR_SIGMA_DEFAULT,
                        help='Trigger threshold in standard deviation.')

    parser.add_argument('-l', '--length',
                        type=int,
                        default=DATA_LEN_DEFAULT,
                        help='Data length.')

    parser.add_argument('-c', '--count',
                        type=int,
                        default=THR_COUNT_DEFAULT,
                        help='Threshold count.')

    parser.add_argument('-r', '--rate',
                        type=int,
                        default=RATE_KSPS_DEFAULT,
                        help='Sampling rate in kHz.')

    parser.add_argument('-p', '--power',
                        type=int,
                        default=1,
                        help='# of ch used for each comm.(<= max_ch in FPGA). default=1')

    parser.add_argument('--position',
                        type=int,
                        default=TRIG_POS_DEFAULT,
                        help='Trigger position.')

    parser.add_argument('-m', '--pre_length',
                        type=int,
                        default=PREMEAN_LEN_DEFAULT,
                        help='Pre-measurement length to estimate data fluctuation.')

    parser.add_argument('-ip', '--ip_address',
                        type=str,
                        default=IP_ADDRESS_DEFAULT,
                        help='IP-v4 address of target SiTCP. (default=192.168.10.16)')

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

    args = parser.parse_args()

    fpga = FPGAControl()
    trig_ch = fpga.trig_ch


    # Argument validation.
    if not trig_ch > 0:
        raise TrgError(f'Firmware version: {fpga.info.version}' +
                       'Trigger function is disabled!')

    if len(args.freqs) > trig_ch:
        raise TrgError('Number of tones exceeds capability of the firmware.' +
                       f'required: {len(args.freq)}/capable: {trig_ch}')

    if len(args.freqs) == 0:
        raise TrgError('No tones specified.')

    dds_f_megahz = args.freqs

    if 200000 % args.rate != 0:
        raise TrgError('Sampling rate (kHz) should divide 200000 evenly.')

    if args.fname is None:
        fname  = 'tod_trg'

        for freq in dds_f_megahz:
            fname += f'_{freq:+08.3f}MHz'

        fname += strftime('_%Y-%m%d-%H%M%S')
        fname += '.rawdata'
    else:
        fname = args.fname

    fname = Path(fname)

    if fname.exists():
        raise TrgError(f'{fname:s} exists.')

    tone_conf = ToneConf(fpga.max_ch, dds_f_megahz,
                         phases=args.phase,
                         amps=args.amplitude,
                         power=args.power)

    measure_trg(fpga,
                tone_conf,
                args.length,
                args.threshold,
                args.count,
                args.rate,
                args.position,
                args.pre_length,
                fname)

if __name__ == '__main__':
    main()
