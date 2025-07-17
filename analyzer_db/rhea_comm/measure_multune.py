import sys
import numpy as np
from time    import strftime
from argparse import ArgumentParser

path_rc = '/home/tomonaga/scripts/analyzer_db'
sys.path.append(path_rc)


import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-a_sta", "--amp_start", help="start of amplitude", default = 0.1, type = float)
parser.add_argument("-a_fin", "--amp_stop", help="stop of amplitude", default = 1.0, type = float)
parser.add_argument("-a_int", "--amp_step", help="step of amplitude", default = 0.1, type = float)
parser.add_argument("-kl", "--kids_list", help="path to kids list", type=str)
parser.add_argument("-w", "--sweep_width", help="full width of sweep range[MHz]", default = 3, type=float)
parser.add_argument("-s", "--step", help="sampling step[MHz]", default = 0.01, type=float)
#parser.add_argument("-f", "__fine_name", help="name of generated file", default=None, type=float)
parser.add_argument("-ip", "--ip_address", help='ip address of fpga', type=str)
parser.add_argument("-port", "--port", help="port of sg", type=str)
args = parser.parse_args()

#path_kl = '/home/tomonaga/scripts/analyzer_db/mkid_pylibs'
from kidslist import KidsList

klist_path = args.kids_list
klist = KidsList(klist_path)
print(klist.kids_freqs)

a_sta = args.amp_start
a_sto = args.amp_stop
a_ste = args.amp_step
w = args.sweep_width
s = args.step
ip = args.ip_address
port = args.port


import os

os.system("python3 sg_manager.py --port " + str(port) + " -p on -f " + str(klist.sg_freq/1e9) + "GHz")

#GB_num, amp
try:
    
    for i in range(len(klist.kids_index)):
            freq = klist.kids_freqs[i]
            p = klist.kids_power[i]
            a = np.arange(a_sta, a_sto, a_ste)
            
            for j, amp in enumerate(a):
                fname  = 'mulswp'
                fname += '_GB_0' + ip[-1]
                fname += f'_{w:+08.3f}MHzWidth'
                fname += f'_{s:+08.3f}MHzStep'
                f_cen = klist.kids_freqs[i]
                fname += f'_{f_cen:+08.3f}MHz'
                #if f_off is not None:
                #    fname += f'_{f_off:+08.3f}MHzOffTone'
                fname += f'_amp_{amp:+08.2f}'
                fname += strftime('_%Y-%m%d-%H%M%S')
                fname += '.rawdata'
                os.system("python3 /home/tomonaga/scripts/analyzer_db/rhea_comm/measure_mulswp.py "\
                    + str(freq)\
                    + " -w " + str(w)\
                    + " -s " + str(s)\
                    + " -f " + fname\
                    + " -p " + str(p)\
                    +" --amplitude " + str(amp)\
                    + " -ip " + str(ip))
except KeyboardInterrupt:
    print("Keyboard interrupt exception caught")