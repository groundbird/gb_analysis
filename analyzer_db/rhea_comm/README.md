rhea_comm
====
Python script to control the rhea data.

## Description

The rhea_comm communicates with the rhea.
Executable files are:
- `measure_swp.py` : taking data with sweeping RHEA frequency
- `measure_mulswp.py` : taking data with sweeping RHEA multiple frequencies
- `measure_sgswp.py` : taking data with sweeping SG frequency
- `measure_tod.py` : taking data with fixing RHEA and SG frequencies with given SPS rate.
- `measure_trg.py` : taking data when firing the trigger with fixing RHEA and SG frequencies.
- `reader_swp.py` : converting taken rawdata file to readable text
- `reader_tod.py` : converting taken rawdata file to readable text
- `packet_reader.py` : functions to read rawdata packets.
- Modules to be used for the analysis are written in `lib_read_rhea.py`.
- for kcu105 control : `adc_, dac_, dds_, debug_, ds_, info_, iq_, raw_, snap_, trg_ setting.py`, `fpga_control.py`, `rbcp .py, _comm.py`, `tcp.py`

## Used Packages

```js
import math
import numpy
import os
import struct
import sys
import time
import socket
import serial
import argparse
import enum
```

If you want to use `measure_sgswp_serial.py`, you should import `sg` module.
It is good to put a symbolic link of the `sg` module directory here.
```bash
$ ln -s ../sg .
```

## Python Version

`python 3.6 or later`

## Usage
before using measure_*.py, you can check the usage by option '-h'.

ex)
```
 > measure_tod.py -h
 >
 > usage: measure_tod.py [-h] [-f FNAME] [-l LENGTH] [-r RATE] [-p POWER]
 >                       [-a AMPLITUDE [AMPLITUDE ...]] [-ip IP_ADDRESS]
 >                       freqs [freqs ...]
 >
 > positional arguments:
 >   freqs                 list of target frequency_MHz.
 >
 > optional arguments:
 >   -h, --help            show this help message and exit
 >   -f FNAME, --fname FNAME
 >                         output filename.
 >                         (default=tod_INPUTFREQ_RATE_DATE.rawdata)
 >   -l LENGTH, --length LENGTH
 >                         data length. (default=10000)
 >   -r RATE, --rate RATE  sample rate of measurement (kSPS). (default=1kSPS)
 >   -p POWER, --power POWER
 >                         # of ch used for each comm.(<= max_ch in FPGA).
 >                         default=-1
 >   -a AMPLITUDE [AMPLITUDE ...], --amplitude AMPLITUDE [AMPLITUDE ...]
 >                         list of amplitude scale (0 to 1.0). # of amplitude
 >                         scale must be same as # of input freqs.
 >   -ip IP_ADDRESS, --ip_address IP_ADDRESS
 >                         IP-v4 address of target SiTCP. (default=192.168.10.16)
```

