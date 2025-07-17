#!/usr/bin/env/ python3

import numpy as np
import sys
# add kidslist path
#sys.path.append("/Users/suenoyoshinori/CMB/gitlab/analyzer_db/")
from kidslist import KidsList
import math


# phase optimize-------
## kitayoshi algorithm
def sigma(val):
    result=0
    for i in range(val):
        result += i + 1
    return result

def phase_opt(num_tone):
    phase_single =[float(sigma(i)*2*math.pi/(num_tone)) for i in range(num_tone)]
    return phase_single
# --------------------------

# cal crest factor-------
def sinf(t, freq, amp, fai):
    return amp*np.sin(2*np.pi*freq*t+fai)

def cosf(t, freq, amp, fai):
    return amp*np.cos(2*np.pi*freq*t+fai)

def cal_cre(klist, phase_list):
    t = np.arange(0,5e-7, 1e-10)
    num_tone = len(phase_list)
    rf_cos_mix = np.zeros(len(t))
    rf_freqlist = np.zeros(num_tone)
    for i,IF_freq in enumerate(klist.tone_freqs):
        rf_freqlist[i] = IF_freq *1e6 + klist.sg_freq
    
#    if len(klist.kids_power) == 1:
    for i,phase in enumerate(phase_list):
        cos = cosf(t, rf_freqlist[i], 1, phase)/num_tone
        rf_cos_mix += cos

## future plan? ------------------
#    elif len(klist.tone_phases) == len(klist.tone_power):
#        for i, (power, phase) in enumerate(zip(power_list, klist.tone_phases)):
#            cos = power * cosf(t, rf_freqlist[i], 1, phase)/num_tone
#            rf_cos_mix += cos
#    else:
#        print('Different between length of phase_list and power_list.')
#---------------------------------
    print('RF frequency lisf')
    print(rf_freqlist)

    max_rf_cos = np.abs(rf_cos_mix).max()
    rms_rf_cos = np.sqrt(np.square(rf_cos_mix).mean())

    rf_cos_crest = max_rf_cos/rms_rf_cos
    return rf_cos_crest
#----------------------------  

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('path',
                        type=str,
                        help='path of kids.list.')

    args = parser.parse_args()
    
    klist = KidsList(args.path)
    num_tone = len(klist.kids_freqs + klist.blinds_freqs)
    
    phase_list = phase_opt(num_tone)

    # make kids_phase and blinds_phase
    klist.kids_phases = [None]*len(klist.kids_freqs)
    klist.blinds_phases = [None]*len(klist.blinds_freqs)
    for i, kid_id in enumerate(klist.kids_index):
        klist.kids_phases[i] = phase_list[kid_id]
    for i, blind_id in enumerate(klist.blinds_index):
        klist.blinds_phases[i] = phase_list[blind_id]
    
    print('phase list')
    print(phase_list)
    print('kids_phase list')
    print(klist.kids_phases)
    print('blinds_phase list')
    print(klist.blinds_phases)

    klist._kldict['kids_phases']= klist.kids_phases
    klist._kldict['blinds_phases']= klist.blinds_phases


    crest_factor = cal_cre(klist,phase_list)
    print('crest factor')
    print(crest_factor)
    klist._kldict['crest_factor']=crest_factor

# Overwrite phase_list and crest factor
    import json
    
    with open(args.path, 'w') as a:
        json.dump(klist._kldict,a,indent=4)
