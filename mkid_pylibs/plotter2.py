import matplotlib.pyplot as plt
import numpy as np
from . import misc


def diff_psd(on, off, kind = 'data'):
    on.calcpsd(kind = kind, dofit = False)
    off.calcpsd(kind = kind, dofit = False)
    plt.figure(figsize = (14,10))
    plt.rcParams['font.size'] = 18
    plt.semilogx(on.psd.f, 10*np.log10(on.psd.amplitude),'.:', c= 'r', label = 'on_amp')
    plt.semilogx(on.psd.f, 10*np.log10(on.psd.phase),'.:', c = 'b', label = 'on_phase')
    plt.semilogx(off.psd.f, 10*np.log10(off.psd.amplitude),'.:', c = 'm', label = 'off_amp')
    plt.semilogx(off.psd.f, 10*np.log10(off.psd.phase),'.:', c= 'c', label = 'off_phase')
    plt.xlabel('freqency [Hz]')
    plt.ylabel('PSD [dBc/Hz]')
    plt.title('different of {} PSD'.format(kind))
    plt.legend(loc = 'best')

