# -*- coding: utf-8 -*-
#### mkid simple peak serach code
"""
MKID peak search code
Main function: find_peaks()

The ampl data points are smoothed by savgol_fitter().
The second derivatives of smoothed ampl data (`deriv2`) is used.
All peaks larger than `threshold` are considered for the fit. (default: 2-sigma of `deriv2`)
Each peak is fitted with lorentzian (breit-wigner) distribution.
The fitted peak height is expected to be : `( height - bkgd ) / bkgd < maxratio`
If fitted peak height is small, refit with wider range will be performed.
The quality-factor is also required to be larger than `minq`.
The fit result with closet resonant frequency to `fc` is returned.

Interface functions:
 - search_peaks(data, fc, Q_search=100, S21min=1, plot=False)

 - search_peak(data, Q_search=100, S21min=1, plot=False):
   No resonant frequency is given.
   The fit result with most centered frequency is returned.

"""

import numpy as np
import scipy
from scipy import signal
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

from misc import MKIDDataException, oddify

def find_peaks(freq, ampl, fc=-1,
               smooth=1, polyorder=0, deriv=0,
               threshold=None, nsigma=None,
               minq=10000, maxratio=0.5):
    """
    Search peaks in data.

    :param freq: frequency data array
    :param ampl: amplitude data array
    :param fc: resonant frequency [GHz]
    :param smooth: float object for smoothing amplitude. defaut = 1 (nosmooth)
    :param polyorder: int object for smoothing polynominal order. default = 0 (nosmooth)
    :param deriv: int object for returning differential order in smoothing.
    :param threshold: float object of deriv2 peak threshold
    :param nsigma: int object for setting `threshold = nsigma * sigma` instead of the threshold value.
    :param minq: required quality-factor as `q > minq`
    :param maxratio: required peak-height as `( height - bkgd ) / bkgd < maxratio`

    :return: list of fitte parameter's dictionary;
             'Q','f0','f0ind','dl','dr','a_off','a_on'
    """

    dfreq = freq[1] - freq[0]

    deriv2 = savgol_filter(ampl, oddify(smooth), polyorder, 2, delta=dfreq)
    ampl_s = savgol_filter(ampl, oddify(smooth), polyorder, 0, delta=dfreq)

    if threshold is None:
        d2sigma   = np.std(deriv2)
        d2sigma   = np.std(deriv2[np.abs(deriv2)<2*d2sigma])
        threshold = nsigma*d2sigma

    # collect all peaks
    kid_indices = []
    for k in np.where(deriv2 > threshold)[0]:
        if k == 0 or k == len(deriv2)-1: continue
        if deriv2[k-1] <= deriv2[k] and deriv2[k] >= deriv2[k+1]:
            kid_indices.append(k)

    if not kid_indices: return []

    # traverse to zero-crossing
    kids = []
    nbadq = 0
    nbadd = 0

    for k in kid_indices:
        l, r = k, k
        while l > 0           and deriv2[l] > 0: l -= 1
        while r < len(deriv2) and deriv2[r] > 0: r += 1
        w = (r - l + 1) * dfreq
        w = w * 6.0 / np.sqrt(3) # convert to FWHM
        l = int((l-k) * 6.0 / np.sqrt(3) + k)
        r = int((r-k) * 6.0 / np.sqrt(3) + k)
        if l < 0: l = 0
        if r >= len(freq): r = len(freq) - 1
        q0 = freq[k] / w
        if q0 < 0.3*minq:
            q0 = minq

        f1, q1, d1, bg = fitLorentzian(freq[l:r], ampl[l:r], freq[k], q0)
        if (bg-d1)/bg > maxratio and l-10 >= 0 and r+10 < len(freq):
            f1, q1, d1, bg = fitLorentzian(freq[l-10:r+10], ampl[l-10:r+10], freq[k], q0)

        if q1 < minq:
            nbadq += 1
            continue
        if (bg-d1)/bg > maxratio:
            nbadd += 1
            continue

        kids.append((f1, q1, d1, bg))

    if nbadq > 0 and len(kid_indices)>1:
        print(f'removed {nbadq} peaks with bad Q')

    if nbadd > 0 and len(kid_indices)>1:
        print(f'removed {nbadd} peaks with bad S21min')

    # sort by frequency
    kids.sort()

    ## pick up a peak which is closest to the carrier frequency (when fc in freq range)
    if len(kids)>0 and fc>=freq[0] and fc<=freq[-1]:
        idx = np.argmin( abs(np.array(kids).T[0]-fc) )
        kids = [kids[idx]]

    for i, (f, q, depth, bg) in enumerate(kids):
        f0ind = np.argmin(abs(freq - f))
        w     = f/q
        dl    = w/2.0/dfreq
        dr    = w/2.0/dfreq
        bg_l  = ampl_s[max(int(f0ind - 3*dl), 0)]
        if int(f0ind-3*dl)<0:
            bg_r = ampl_s[0]
        else:
            bg_r  = ampl_s[min(int(f0ind - 3*dl), len(freq) - 1)]
        a_off = (bg_l + bg_r)/2.0
        a_on  = ampl_s[f0ind]
        kids[i] = {'Q': q, 'f0': f, 'f0ind': f0ind,
                  'dl': int(dl), 'dr': int(dr),
                   'a_off': a_off, 'a_on': a_on}

    return kids

def search_peaks(data, **kws):
    return find_peaks(data.x, data.amplitude, **kws)

def search_peak_center(data, **kws):
    peaks = find_peaks(data.x, data.amplitude, **kws)
    if not peaks:
        raise MKIDDataException('peak find failure')

    center  = (data.x[0] + data.x[-1])/2.0
    minind  = 0
    mindist = abs(peaks[0]['f0'] - center)

    for i in range(1, len(peaks)):
        if mindist > abs(peaks[i]['f0'] - center):
            minind  = i
            mindist = abs(peaks[i]['f0'] - center)

    return peaks[minind]

def fitLorentzian(freq, ampl, f0, q0):
    """
    Fit data with lorenzian curve.

    :param freq: a 1-D array of frequency
    :param ampl: a 1-D array of amplitude
    :param f0: initial parameter for center frequency
    :param q0: initial parameter for quality factor

    :return: (fc, q, d, bg)

    - **fc**: fit of center frequency
    - **q**: fit of quality factor
    - **d**: fit of amplitude for Lorentzian curve
    - **bg**: fit of constant background level
    """
    def f(x):
        (a, b, c, d) = x # background, amplitude, 2/FWHM, freq. center
        y = a + b / (((freq-d)*c)**2 + 1)
        return y - ampl
    def fprime(x):
        (a, b, c, d) = x
        g = np.zeros((len(freq)), 4)
        g[:, 0] = 1.0
        g[:, 1] = 1.0 / (1.0 + ((freq - d)*c)**2)
        g[:, 2] = -2.0 * b * c * (freq - d)**2 / (1.0 + ((freq - d)*c)**2)**2
        g[:, 3] =  2.0 * b * c**2 * (freq - d) / (1.0 + ((freq - d)*c)**2)**2

    a = np.median(ampl)
    b = -0.8 * a
    c = 2.0 * q0 / f0
    d = f0
    x0 = np.array([a, b, c, d])
    x1 = scipy.optimize.leastsq(f, x0)[0]
    (a, b, c, d) = x1
    f = d
    q = abs(c * d / 2.0)

    return (f, q, -b, a)

def analytical_search_peak(data):

    y0 = data.amplitude
    x0 = data.x

    # linear bkgd
    a0 = (y0[-1] - y0[0])/(x0[-1] - x0[0])
    b0 = y0[0] - a0*x0[0]

    f0ind = np.argmin(y0 - (a0*x0+b0))
    f = x0[f0ind]
    flind = np.argmin(np.abs(y0[f0ind]/2.0 + (a0*x0[f0ind]+b0)/2.0 - y0[:f0ind]))
    frind = np.argmin(np.abs(y0[f0ind]/2.0 + (a0*x0[f0ind]+b0)/2.0 - y0[f0ind:])) + f0ind
    dl = f0ind - flind
    bgl = y0[np.amax([int(f0ind-3*dl),0])]
    dr = frind - f0ind
    bgr = y0[np.amin([int(f0ind+3*dr),len(y0)-1])]
    a_on = y0[f0ind]
    a_off = (bgl+bgr)/2.0
    q = f/(x0[frind]-x0[flind])

    return {'Q': q, 'f0': f, 'f0ind': f0ind,
            'dl': int(dl), 'dr': int(dr),
            'a_off': a_off, 'a_on': a_on}
    
def analytical_search_peak2(data, dep = 3):

    y0 = data.amplitude
    x0 = data.x

    # linear bkgd
    a0 = (y0[-1] - y0[0])/(x0[-1] - x0[0])
    b0 = y0[0] - a0*x0[0]
    
    #flag = (y0 - np.mean(y0)) < 0
    #print(y0[np.argmin(y0)])
    #print(- y0[np.argmin(y0)] - (y0[np.argmax(y0)] - y0[np.argmin(y0)])/5)
    flag = (y0 - y0[np.argmin(y0)] - (y0[np.argmax(y0)] - y0[np.argmin(y0)])/dep) < 0
    #print(x0[flag])
    print(np.diff(x0[flag]))
    max_ind = np.argmax(np.diff(x0[flag]))
    print(f'max_ind : {max_ind}')
    #print(x0[flag][max_ind])
    f1ind = np.argmin(y0[flag][:max_ind+1])
    f2ind = np.argmin(y0[flag][max_ind+1:-1])
    f1 = x0[flag][:max_ind+1][f1ind]
    f2 = x0[flag][max_ind+1:-1][f2ind]
    print(f'f1 : {f1}')
    print(f'f2 : {f2}')
    
    # first KID
    flind1 = np.argmin(np.abs(y0[f1ind]/2.0 + (a0*x0[f1ind]+b0)/2.0 - y0[:f1ind]))
    frind1 = np.argmin(np.abs(y0[f1ind]/2.0 + (a0*x0[f1ind]+b0)/2.0 - y0[f1ind:])) + f1ind
    dl1 = f1ind - flind1
    dr1 = frind1 - f1ind
    q1 = f1/(x0[frind1]-x0[flind1])
    #print('first KID')
    #print(dl1)
    #print(dr1)
    #print(q1)
    
    
    # second KID
    flind2 = np.argmin(np.abs(y0[f2ind]/2.0 + (a0*x0[f2ind]+b0)/2.0 - y0[:f2ind]))
    frind2 = np.argmin(np.abs(y0[f2ind]/2.0 + (a0*x0[f2ind]+b0)/2.0 - y0[f2ind:])) + f2ind
    dl2 = f2ind - flind2
    dr2 = frind2 - f2ind
    q2 = f2/(x0[frind2]-x0[flind2])
    #print('second KID')
    #print(dl2)
    #print(dr2)
    #print(q2)
    
    
    # on/off amplitude values
    a_on1 = y0[f1ind]
    a_on2 = y0[f2ind] 
    bgl = y0[np.amax([int(f1ind-3*dl1),0])]
    bgr = y0[np.amin([int(f2ind+3*dr2),len(y0)-1])]       
    a_off = (bgl+bgr)/2.0

    return {'Q1': q1, 'Q2': q2, 'f1': f1, 'f2': f2, 'f1ind': f1ind, 'f2ind': f2ind,
            'a_off': a_off, 'a_on1': a_on1, 'a_on2': a_on2}


