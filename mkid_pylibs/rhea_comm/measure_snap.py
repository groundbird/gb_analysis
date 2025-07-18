#!/usr/bin/env python3


## config
cnt_step   = 100000   # byte
cnt_finish = 10000000 # It is larger than FPGA FIFO size.
source_list = ['DDS', 'DAC', 'ADC', 'IQ', 'TEST', 'PHASE']
source_default = 'IQ'


## constant
from fpga_control import FPGAControl
fpga = FPGAControl()
MAX_CH = fpga.max_ch

## check firmware
from sys import stderr
if not fpga.en_snap:
    print(f'Firmware version: {fpga.info.version:d}', file=stderr)
    print('Snapshot function is disable!', file=stderr)
    exit(1)
    pass

## main
def measure_snap(dds_f_MHz, source, channel, fname):
    from math    import floor
    from struct  import pack

    print('SNAPSHOT')
    for i, freq in enumerate(dds_f_MHz):
        print(f'ch{i:03d}: {freq} MHz')
        pass

    fpga.init()
    fpga.dds_setting.set_freqs([freq * 1e6 for freq in dds_f_MHz])
    fpga.snap_setting.set_src(source, channel)

    f = open(fname, 'wb')

    ## write header
    for ch, freq in enumerate(dds_f_MHz):
        dummy_packet  = b'\xff'
        dummy_packet += b'\xff' * 4
        dummy_packet += pack('b', -1-ch)
        freq_Hz = int(floor(freq * 1e6 + 0.5))
        freq_packet   = pack('>i', freq_Hz)
        dummy_packet += freq_packet * 2
        dummy_packet += b'\xee'
        f.write(dummy_packet)
        pass

    cnt = 0
    cnt_print = cnt_step
    empty_count = 0

    fpga.tcp.clear()
    fpga.snap_setting.snap_on()

    while True:
        buff = fpga.tcp.read(1024)
        f.write(buff)
        if not len(buff): break
        cnt += len(buff)
        while cnt >= cnt_print:
            print(cnt_print)
            cnt_print += cnt_step
            pass
        if cnt >= cnt_finish: break
        pass

    f.close()
    fpga.dac_setting.txenable_off()
    print(f'write raw data to {fname}')


if __name__ == '__main__':
    ## check argv
    from sys import argv
    from os.path import isfile
    from time    import strftime

    try:
        args = argv[1:]
        arg_val = []
        fname = None
        source  = source_default
        channel = 0
        while args:
            if args[0] == '-f':
                fname = args[1]
                args = args[2:]
            elif args[0] == '-s':
                source = args[1]
                args = args[2:]
            elif args[0] == '-c':
                channel = int(args[1])
                args = args[2:]
            else:
                arg_val += [args[0]]
                args = args[1:]
                pass
            pass
        if len(arg_val) > MAX_CH:
            raise Exception('exceeds max ch.')
        if len(arg_val) == 0:
            raise Exception('No DDS freq specified.')
        dds_f_MHz = [float(v) for v in arg_val]
        while len(dds_f_MHz) < MAX_CH: dds_f_MHz.append(0.)
        if not source in source_list:
            raise Exception('source error')

        if fname == None:
            fname  = 'snap'
            fname += f'_{source:s}'
            fname += f'_ch{channel:03d}'
            fname += strftime('_%Y-%m%d-%H%M%S')
            fname += '.rawdata'
        if isfile(fname):
            raise Exception(f'{fname:s} exists')
    except:
        print( 'Usage: measure_snap.py [freq0_MHz] [freq1_MHz] ...', file=stderr)
        print( '       option: -f [filename]     (mandantory)', file=stderr)
        print( '                  output filename', file=stderr)
        print(f'               -s [(char)source] (default: {source_default:s})', file=stderr)
        print( '                  ', end=' ', file=stderr)
        for s in source_list[:-1]: print(f'{s:s},', end=' ', file=stderr)
        print(f'{source_list[-1]:s}', file=stderr)
        print( '               -c [(int)channel] (default: 0)', file=stderr)
        print( '                  channel number to be snapped')
        exit(1)

    measure_snap(dds_f_MHz, source, channel, fname)

