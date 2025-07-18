#!/usr/bin/env python3

from math import floor
from time   import sleep
import numpy as np

from packet_reader import read_iq_packet
from fpga_control  import FPGAControl

## constant
RHEACLK_RATE = 200e6 # Hz

RHEADAC_VPP = 1 # V
RHEADAC_DBM = 10*np.log10((((RHEADAC_VPP/2/np.sqrt(2))**2)/50)/1e-3) # dBm

RHEAADC_BIT = 28

def get_packet_size(ntones):
    return 7 + 7 * 2 * ntones

## main
def check_fpga(fpga, max_ch, dds_f_MHz, rate_kSPS=1, power=1, amps=None, phases=None, is_norm=False, verbose=False):

    def vprint(*args,**kwargs):
        if verbose:
            print(*args,**kwargs)
            pass
        return

    buff = None
    try:
        vprint(f'FPGA CHECK')
        sleep(1)

        vprint('-- init()')
        fpga.init()
        if len(dds_f_MHz) == 0: return

        vprint('-- iq.set_read_width()')
        fpga.iq.set_read_width(len(dds_f_MHz))
        vprint('-- tcp.clear()')
        fpga.tcp.clear()

        dds_frq = [freq * 1e6 for freq in dds_f_MHz] * int(power) + [0.] * (max_ch - int(power)*len(dds_f_MHz))
        vprint(f'-- dds.set_freqs()    :  |'+'|'.join([f' {x:+4.3f} ' for x in dds_f_MHz])+'|')
        fpga.dds.set_freqs(dds_frq)

        amps = [1.]*len(dds_f_MHz) if amps is None else amps
        dds_amp = amps * int(power) + [0.] * (max_ch - int(power)*len(amps))
        vprint(f'-- dds.set_amps()     :  |'+'|'.join([f'  {x:1.4f} ' for x in amps])+'|')
        fpga.dds.set_amps(dds_amp)

        phases = [0.]*len(dds_f_MHz) if phases is None else phases
        dds_phs = phases * int(power) + [0.] * (max_ch - int(power)*len(phases))
        vprint(f'-- dds.set_phases()   :  |'+'|'.join([f' {x:+1.4f} ' for x in phases])+'|')
        fpga.dds.set_phases(dds_phs)

        vprint(f'-- dds.set_rate()')
        fpga.ds.set_rate(floor(200.e3 / rate_kSPS + 0.5)) # 1e3 SPS

        vprint('-- iq.time_reset()')
        fpga.iq.time_reset()
        vprint('-- iq.iq_on()')
        fpga.iq.iq_on()

        if len(dds_f_MHz)>0:
            packet_size = get_packet_size(len(dds_f_MHz))
            buff = fpga.tcp.read(packet_size)

    except KeyboardInterrupt:
        print('stop by keyboard input')
    finally:
        vprint('-- iq.iq_off()')
        fpga.iq.iq_off()
        vprint('-- tcp.clear()')
        fpga.tcp.clear()
        vprint('-- dac.txenable_off()')
        fpga.dac.txenable_off()

        if buff is not None:
            vprint('='*80)
            time, data, nrot, sync = read_iq_packet(buff)
            print(f'time = {time}, nrot = {nrot}, sync_off = {sync}, #_IQ_column = {len(data)}, rate = {rate_kSPS:d}kSPS, tone_power = {power:d}/{max_ch:d} = {RHEADAC_DBM + 20*np.log10(power/max_ch):.3}dBm')
            print('FREQ  : '+'  '.join([f'{x:+18.3f}' for x in dds_f_MHz]))
            fact = (2**RHEAADC_BIT) * RHEACLK_RATE / (rate_kSPS*1e3) if is_norm else 1.
            print('S21amp: '+'  '.join([f'{abs(data[2*i] + 1j*data[2*i+1])/fact:1.12e}' for i in np.arange(0,int(len(data)/2))]))
            vprint('='*80)
            pass

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('-ip', '--ip_address',
                        type=str,
                        default='192.168.10.16',
                        help='IP-v4 address of target SiTCP. (default=192.168.10.16)')

    parser.add_argument('freqs',
                        type=float,
                        nargs='*',
                        default = [],
                        help='list of target frequency_MHz.')

    parser.add_argument('-r', '--rate',
                        type=int,
                        default=1,
                        help=f'sample rate of measurement (kSPS). (default=1kSPS)')

    parser.add_argument('-p', '--power',
                        type=int,
                        default=1,
                        help=f'# of ch used for each comm.(<= max_ch in FPGA). default=1')

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

    parser.add_argument('-n', '--normalize',
                        action='store_true',
                        help='dump normalized S21')

    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Verbose mode.')

    args = parser.parse_args()
    ip          = args.ip_address
    freqs       = args.freqs
    rate_kSPS   = args.rate
    power       = args.power
    amps        = args.amplitude
    phases      = args.phase
    norm        = args.normalize

    import subprocess

    try:
        status,result = subprocess.getstatusoutput('ping -c1 '+ip)
        print('='*80)
        print(result)
        print('='*80)
        if status == 1:
            print(f"No ping resoponse from {ip}...")
            raise TimeoutError
        fpga = FPGAControl(ip_address=ip)
        max_ch = fpga.max_ch
        dds_f_MHz = []
        if len(freqs)>0:
            if len(freqs) > max_ch:
                raise Exception(f'too many frequencies! = {len(freqs)} / {max_ch}')
            if power>0 and len(freqs)*power > max_ch:
                raise Exception(f'exceeding max # of channels = {len(freqs)}*{power} > {max_ch}')

            def two_div(val):
                if val < 2: return 0;
                return len(bin(val-1)) -2

            input_len = 2**two_div(len(freqs))
            if power >0 and power*input_len > max_ch:
                raise Exception(f'exceeding max # of channels = {input_len}*{power} > {max_ch}')
            power = floor(float(max_ch) / float(len(freqs))+0.1) if power<0 else power

            if int(RHEACLK_RATE*1e-3) % rate_kSPS != 0:
                raise Exception(f'sampling rate [kSPS] should be a divisor of 2000000: input rate = {rate_kSPS} kSPS')
            dds_f_MHz = [float(v) for v in freqs]
            if amps is not None and len(amps) != len(freqs):
                raise Exception(f'# of input amp must be same as # of input freqs.')
            if phases is not None and len(phases) != len(freqs):
                raise Exception(f'# of input phase must be same as # of input freqs.')
    except TimeoutError:
        print('connection to FAGA failed.')
        print(ip, 'is invalid ip address.')
        exit(1)
        pass
    except Exception as e:
        print(e)
        exit(1)
        pass

    check_fpga(fpga        = fpga,
               max_ch      = max_ch,
               dds_f_MHz   = dds_f_MHz,
               rate_kSPS   = rate_kSPS,
               power       = power,
               amps        = amps,
               phases      = phases,
               is_norm     = norm,
               verbose     = args.verbose)


