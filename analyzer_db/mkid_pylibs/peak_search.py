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

from .misc import MKIDDataException, oddify

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
'''
def analytical_search_peak(data):

    y0 = data.amplitude
    ytmp = scipy.signal.medfilt(y0,3)
    spike_index = np.where(np.abs((ytmp - y0)/ytmp)>0.1)[0]
    for i in spike_index:
        if i == 0 or i == len(y0)-1: continue
        y0[i] = (y0[i-1]+y0[i+1])/2

    x0 = data.x

    # linear bkgd
    a0 = (y0[-1] - y0[0])/(x0[-1] - x0[0])
    b0 = y0[0] - a0*x0[0]

    def get_peakindex(f,amp,std):
        df = f[1]-f[0]
        width = (int(50*1e3/df+0.1),int(5*1e6/df+0.1)) # 50k -- 5MHz
        prominence = (10**(int(np.log10(std) - 1))/2, None)
        distance = int(50*1e3/df+0.1)
        ret = scipy.signal.find_peaks(-amp,width=width,prominence=prominence,distance=distance)[0]
        if len(ret)==0:
            amp2 = scipy.signal.medfilt(amp,3)
            ret = scipy.signal.find_peaks(-amp2,width=width,prominence=prominence,distance=distance)[0]
        return ret

    y0subt = y0 - (a0*x0+b0)
    peak_inds = get_peakindex(x0,y0subt,np.nanstd(y0subt))
#    tmparr = x0.mean()-x0[peak_inds]
#    f0ind = peak_inds[np.argmin([np.abs(a)*(1.5 if a>0 else 1) for a in tmparr])]
    f0ind = peak_inds[np.argmin(np.abs(x0.mean()-x0[peak_inds]))]
#    f0ind = peak_inds[np.argmin(y0subt[peak_inds]-y0subt.mean())]


    #f0ind = np.argmin(y0 - (a0*x0+b0))
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
'''
# used analytical_search_peak below because scipy version is old.
def analytical_search_peak(data):

    y0 = data.amplitude
    x0 = data.x

    # linear bkgd
    a0 = (y0[-1] - y0[0])/(x0[-1] - x0[0])
    b0 = y0[0] - a0*x0[0]

    f0ind = np.argmin(y0 - (a0*x0+b0))
    f = x0[f0ind]
    flind = np.argmin(np.abs(y0[f0ind]/2.0 + (a0*x0[f0ind]+b0)/2.0\
 - y0[:f0ind]))
    frind = np.argmin(np.abs(y0[f0ind]/2.0 + (a0*x0[f0ind]+b0)/2.0\
 - y0[f0ind:])) + f0ind
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
