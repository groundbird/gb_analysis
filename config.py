
"""
Configuration file for KID data processing system.

This module contains all chip-specific configurations and system constants
used for KID resonance fitting and data processing.
"""

import numpy as np

# System Configuration
CHIPS = np.array(['3A', '2A', '3B', '1A', '2B', '1B', '220'])
DAQS = np.array(['GB01', 'GB02', 'GB03', 'GB04', 'GB05', 'GB06', 'GB07'])
SAVEDIR = './'



def get_fitconfig(chip):
    # configuration for fitting each detectors
    # See mkid_pylibs_sueno for each configuration.
    if True:
        if chip == '1A':
            rangeind = [15]
            freqranges = [(4.886e9, 4.8885e9)]
            twokidind = [5, 6]
            twokidfitter = ['gaolinbg2f', 'gaolinbg2l']
            initind = []
            fitinit = []
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '1B':
            rangeind = [10, 17, 19, 20, 21, 22] # KID2,3,4 are degenerate.
            freqranges = [(4.8105e9, 4.8125e9), (4.878e9, 4.881e9), (4.888e9, 4.890e9), (4.890e9, 4.892e9), (4.894e9, 4.8964e9), (4.8965e9, 4.899e9)]
            twokidind = [9, 17, 18]
            twokidfitter = ['gaolinbg2f', 'gaolinbg2f', 'gaolinbg2l']
            initind = []
            fitinit = []
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '2A':
            rangeind = []
            freqranges = []
            twokidind = []
            twokidfitter = []
            initind = []
            fitinit = []
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '2B':
            rangeind = []
            freqranges = []
            twokidind = []
            twokidfitter = []
            initind = []
            fitinit = []
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '3A':
            rangeind = [7, 8, 9, 10, 13, 14]
            freqranges = [(4.736e9, 4.7385e9), (4.736e9, 4.7385e9), (4.7385e9, 4.741e9), (4.7385e9, 4.741e9), (4.7825e9, 4.785e9), (4.7825e9, 4.785e9)]
            #twokidind = [7, 8, 9, 10, 13, 14]
            #twokidfitter = ['gaolinbg2f', 'gaolinbg2l', 'gaolinbg2f', 'gaolinbg2l', 'gaolinbg2f', 'gaolinbg2l']
            twokidind = [7, 8, 10, 13, 14]
            twokidfitter = ['gaolinbg2f', 'gaolinbg2l', 'gaolinbg2l', 'gaolinbg2f', 'gaolinbg2l']
            #initind = [7, 8, 9, 10, 13, 14]
            initind = [7, 8, 10, 13, 14]
            fitinit = [{'fr1' : 4.7373e9, 'fr2' : 4.7376e9, 'Qr1': 15000, 'Qr2': 30000, 'Qc1': 15000, 'Qc2': 30000, 'arga': 0, 'absa': 0.01, 'tau': 0, 'phi01': 0, 'phi02': 0, 'c': 0 }, # 7
                       {'fr1' : 4.7373e9, 'fr2' : 4.7376e9, 'Qr1': 15000, 'Qr2': 30000, 'Qc1': 15000, 'Qc2': 30000, 'arga': 0, 'absa': 0.01, 'tau': 0, 'phi01': 0, 'phi02': 0, 'c': 0 }, # 8
                       #{'fr1' : 4.7394e9, 'fr2' : 4.7399e9, 'Qr1': 10000, 'Qr2': 30000, 'Qc1': 17000, 'Qc2': 40000, 'arga': 0, 'absa': 0.01, 'tau': 0, 'phi01': 0, 'phi02': 0, 'c': 0 }, # 9
                       {'fr1' : 4.7395e9, 'fr2' : 4.7399e9,'Qr1': 16000, 'Qr2': 40000, 'Qc1': 17000, 'Qc2': 40000}, # 10
                       {'fr1' : 4.7836e9, 'fr2':4.7842e9, 'Qr1': 15000, 'Qr2': 35000, 'Qc1': 10000, 'Qc2': 20000}, # 13
                       {'fr1' : 4.7836e9, 'fr2':4.7842e9, 'Qr1': 15000, 'Qr2': 35000, 'Qc1': 10000, 'Qc2': 20000}] # 14
            depind = [8, 9]
            dep = [2, 2]
            #skipind = [7, 8, 9]
            #guessskip = [True, True, True]
            skipind = [7]
            guessskip = [True]
            
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '3B':
            #rangeind = [0, 12]
            #freqranges = [(4.6858e9, 4.6868e9), (4.819e9, 4.821e9)]
            rangeind = [0]
            freqranges = [(4.68575e9, 4.6868e9)]

            twokidind = []
            twokidfitter = []
            initind = [0]
            fitinit = [{'fr' : 4.6862e9, 'Qr' : 15000, 'Qc' : 6000, 'phi0' : 0.1}]
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        elif chip == '220':
            rangeind = [0, 2, 15] 
            freqranges = [(4.854e9, 4.857e9), (4.867e9, 4.869e9), (4.993e9, 4.995e9)]
            twokidind = [14]
            twokidfitter = ['gaolinbg2f']
            initind = []
            fitinit = []
            depind = []
            dep = []
            skipind = []
            guessskip = []
            return rangeind, freqranges, twokidind, twokidfitter, initind, fitinit, depind, dep, skipind, guessskip
        else:
            print('Error : chip:' + chip + 'is not matched.')