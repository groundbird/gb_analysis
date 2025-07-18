#!/usr/bin/env python3


## constant
from fpga_control  import fpga_control
fpga = fpga_control()
MAX_CH = fpga.MAX_CH
packet_size = 7 + 7 * 2

from sg import sg_manager

## main
def measure_sgswp(f_start, f_end, f_step, fname, power):
    from packet_reader import read_iq_packet
    from numpy  import arange, floor
    from struct import pack
    from time   import sleep

    print('SWEEP MEASUREMENT')
    print(f'SwpPower: {power:d}/{MAX_CH:d}')
    sleep(1)

    fpga.init()
    fpga.iq.set_read_width(1)
    cnt_finish = packet_size * 10

    # SG initialization
    quicksyn = sg_manager.QuickSyn()
    sleep(1.0)

    f = open(fname, 'wb')
    fpga.tcp.clear()

    for freq in arange(f_start, f_end, f_step):
        print(f'{freq:8.3f} MHz')
        freq_kHz = int(floor(freq * 1e3 + 0.5))
        freq_str = f'{freq+50}MHz'
        quicksyn.set_freq_str(freq_str)

        sleep(0.1)
        cnt = 0

        #sleep(1)
        while True:
            fpga.dds.set_freqs([-50.0*1e6] * power + [0.] * (MAX_CH - power))
            fpga.iq.time_reset()
            fpga.iq.iq_on()

            try:
                while True:
                    time, data = read_iq_packet(fpga.tcp.read(packet_size))
                    if time == 0: break
                    pass
                break
            except:
                fpga.iq.iq_off()
                fpga.tcp.clear()
                print('... retry read_iq_packet()')
                continue
            pass

        freq_pack  = (b'\x00' * 3) if freq_kHz >= 0 else (b'\xff' * 3)
        freq_pack += pack(">i", freq_kHz)
        dummy_packet  = b'\xff'        # header
        dummy_packet += b'\x00' * 5    # time
        dummy_packet += freq_pack * 2  # freq
        dummy_packet += b'\xee'        # footer
        f.write(dummy_packet)

        while True:
            buff = fpga.tcp.read(min(1024, cnt_finish - cnt))
            f.write(buff)
            if not len(buff): break
            cnt += len(buff)
            if cnt >= cnt_finish: break
            pass

        fpga.iq.iq_off()
        pass

    f.close()
    fpga.dac.txenable_off()
    print(f'write raw data to {fname}')

if __name__ == '__main__':
    ## get argv
    from sys     import argv, stderr
    from os.path import isfile
    from time    import strftime

    try:
        args = argv[1:]
        arg_val = []
        fname = None
        power = 1
        while args:
            if args[0] == '-f':
                fname = args[1]
                args = args[2:]
            elif args[0] == '-p':
                power = int(args[1])
                args = args[2:]
            elif args[0] == '-h':
                print('PRINT HELP::')
                raise
            else:
                arg_val += [args[0]]
                args = args[1:]
                pass
            pass
        if power < 1 or power > MAX_CH:
            print(f'Invalid SwpPower: {power:d}/{MAX_CH:d}', file=stderr)
            raise
        f_start = float(arg_val[0])
        f_end   = float(arg_val[1])
        f_step  = float(arg_val[2])
        if fname is None:
            fname  = 'swp'
            fname += f'_{f_start:+08.3f}MHz'
            fname += f'_{f_end:+08.3f}MHz'
            fname += f'_{f_step:+08.3f}MHz'
            fname += strftime('_%Y-%m%d-%H%M%S')
            fname += '.rawdata'
        if isfile(fname):
            print(f'{fname!r} is existed.', file=stderr)
            raise
    except:
        print( 'Usage: measure_swp.py [freq_start_MHz] [freq_end_MHz] [freq_step_MHz]', file=stderr)
        print( '       option: -f [filename] (mandantory)', file=stderr)
        print( '                  output filename', file=stderr)
        print(f'             : -p [num]      (default: {power})', file=stderr)
        print(f'                  amplitude (<={MAX_CH})', file=stderr)
        exit(1)

    measure_sgswp(f_start, f_end, f_step, fname, power)

