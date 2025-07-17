# -*- coding: utf-8 -*-
#### mkid fit interface
"""
MKID fitting interface

Main function:
fit_onepeak(data, fc, errors=None, nfwhm=3, fitter='gaolinbg', Q_search=10000, plot=False)

"""
import numpy as np
import warnings
import json
import pickle as pickle
import collections

import matplotlib.pyplot as plt

from . import fit

from .misc import MKIDDataException
from . import fitters
#from .peak_search import search_peaks
from .peak_search import analytical_search_peak
import scipy
from scipy import signal

class KidFitResult(object):
    def __init__(self, fitresult=None, fitrange=slice(None)):
        if fitresult is None:
            self._result = fit.fitresult.FitResult()
        else:
            self._result = fitresult
        self.bgfunc      = None
        self.rewindfunc  = None
        self.fitrange    = fitrange
        #self.info        = dict()
    def fitted(self, x):
        return self.eval(x)
    def add_functions(self, funcdict):
        if 'bgfunc' in funcdict:
            self.bgfunc = funcdict['bgfunc']
        if 'rewindfunc' in funcdict:
            self.rewindfunc = funcdict['rewindfunc']
    def rewind(self, x, iq):
        return self.rewindfunc(x, iq, **self.values())
    def bg(self, x):
        return self.bgfunc(x, **self.values())

    @property
    def fitparamdict(self):
        return self.params.valuesdict()

    def dumps(self, header='', **kws):
        from pprint import pformat
        out = dict()
        out['result']     = json.loads(self._result.dumps())
        out['rewindfunc'] = pickle.dumps(self.rewindfunc)
        out['bgfunc']     = pickle.dumps(self.bgfunc)
        out['fitrange']   = pickle.dumps(self.fitrange)
        #out['info']       = pickle.dumps(self.info)

        s = header + '\n\n'
        s += pformat(self.info, **kws) + '\n\n'
        s += self._result.report_str()

        buf = ''
        for l in s.split('\n'):
            buf += ('# ' + l + '\n')
        buf += pformat(out, **kws)
        return buf
        # return json.dumps(out, **kws)
    @classmethod
    def loads(cls, s, **kws):
        r = KidFitResult()
        # dic = json.loads(s, **kws)
        ns = globals()
        ns['inf'] = np.inf
        ns['nan'] = np.nan
        dic = eval(s, ns)
        for k, v in list(dic.items()):
            if k == 'result':
                r._result = r._result.loads(json.dumps(v))
            elif k == 'rewindfunc':
                r.rewindfunc = pickle.loads(v)
            elif k == 'bgfunc':
                r.bgfunc = pickle.loads(v)
            elif k == 'fitrange':
                r.fitrange = pickle.loads(v)
            elif k == 'info':
                r.info = pickle.loads(v)
            else:
                raise RuntimeError('unknown key: %s' % k)
        return r

    def dump(self, fp, **kws):
        return fp.write(self.dumps(**kws))
    @classmethod
    def load(cls, fp, **kws):
        return cls.loads(fp.read(), **kws)

    ## delegate all other methods
    def __getattr__(self, name):
        if name == '_result':
            raise AttributeError()
        return getattr(self._result, name)
    def __getstate__(self):
        return self.dumps()
    def __setstate__(self, s):
        self.__dict__ = self.loads(s).__dict__

def fit_complex(func, fs, iqs, erriqs, params, silent=False):
    """
    Fit a complex-valued function to real-imag plane sweep.

    :param func: a complex-valued function of form f(x, param1, param2, ..., paramN)
    :param fs: a 1-D array of float, of frequency
    :param iqs: a 1-D array of complex data
    :param erriqs: a 1-D array of complex data error
    :param params: a lmfit.Parameters object to fit

    :return: a lmfit.FitResult object
    """

    def convert(x): #complex_to_cartesian2darray
        if x is None: return None
        x = np.atleast_1d(x)
        shape = x.shape
        return np.concatenate([np.real(x), np.imag(x)], axis=len(shape)-1)

    return fit.fit(func, fs, convert(iqs), convert(erriqs), params, convert=convert,silent=silent)

def _fitIQ(data, fitrange, dataerror, params, func, silent=False):
    """
    Fit IQ from SweepData class using given data, fitrange, dataerror, params.

    :param data: a SweepData to fit
    :param fitrange: a slice or tuple, to specify fit range in `data`
    :param dataerror: two-element tuple (Ierr, Qerr) corresponding to error of I and Q, respectively
    :param params: a lmfit.Parameters object to describe parameter
    :param func: a complex-valued function of form f(x, param1, param2, ..., paramN)

    :return: a FitResult object.
    """

    ### fitrange:: cast `s` to slice if it's tuple
    if isinstance(fitrange, slice):
        s = fitrange
    else:
        if isinstance(fitrange[0], np.ndarray): # xrange
            s = slice(fitrange[0][0], fitrange[0][-1]+1)
        else:
            s = slice(*fitrange)

    x = data.x[s]
    iq = data.iq[s]

    if True:
        y = scipy.signal.medfilt(np.abs(iq),3)
        spike_index = np.where(np.abs((y -np.abs(iq))/y)>0.1)[0]
        for i in spike_index:
            if i == 0 or i == len(iq)-1: continue
            #print(f"spike correction: index = {i:4d}, {np.abs(iq[i])} --> {np.abs((iq[i-1]+iq[i+1])/2)} (btw {np.abs(iq[i-1])} and {np.abs(iq[i+1])})")
            iq[i] = (iq[i-1]+iq[i+1])/2

    if dataerror is None:
        erriq = None
    else:
        err_I = dataerror[0]
        if hasattr(err_I,'__iter__') and len(err_I) >1: err_I = err_I[s]
        err_Q = dataerror[1]
        if hasattr(err_Q,'__iter__') and len(err_Q) >1: err_Q = err_Q[s]
        erriq = np.broadcast_to(err_I, iq.shape) + 1j*np.broadcast_to(err_Q, iq.shape)

    fitresult = fit_complex(func, x, iq, erriq, params,silent=silent)
    kidfitresult = KidFitResult(fitresult, s)
    return kidfitresult

#### utility functions
def fitIQ_onepeak(data, errors=None, nfwhm=3, frqrange=[None,None], fitter='gaolinbg', silent=False):
    """
    Fit data with fitter.

    :param data: sweep data to fit
    :param error: error in data (data format depends on fitter)
    :param nfwhm: fit range [FWHM]. if negative, use all data as fit range
    :param fitter: fitter to eb used.
    :param frqrange: list of starting and ending frequencies to be considered

    **fitter** supports:

    - **gao**: a complex function from Jiansong Gao's D thesis (default)
    - **gaolinbg**: a complex function from Jiansong Gao's D thesis
    - **gaolinbg2**: a complex function from Jiansong Gao's D thesis with prefitting
    - **mazinrev**: a complex function from Ben Mazin's D thesis
    - **blank**: a complex function only with cable delay effect or a Fitter object or a tuple (refer fitters.py)
    """
    ## get fitting information
    if type(fitter) == str:
        if fitter not in fitters.all_fitters:
            raise RuntimeError('if you specify fitter by string, it must be one of %s!' % fitter.all_fitters)
        fitter = getattr(fitters, 'fitter_' + fitter)

    func, guess, names, others = fitter
    params = guess(data)

    ## get fitrange
    s0 = None
    s1 = None
    if frqrange[0] is not None:
        s0 = np.where(data.x>frqrange[0])[0][0]
    if frqrange[1] is not None:
        s1 = np.where(data.x<frqrange[1])[0][-1]

    s = slice(s0,s1)
    if nfwhm > 0:
        pdict = analytical_search_peak(data) # assuming one peak in the sweep data
        s = adjust_fitrange(nfwhm, len(data.x), len(names), pdict)

    ## do prefitting and fix some parameters
    if 'prefit_and_fix' in fitter[-1]:
        r_prefit = _fitIQ(data, slice(None), errors, params, func, silent=silent)
        params   = r_prefit.params
        for k in fitter[-1]['prefit_and_fix']:
            params[k].vary = False

    ## add additional exprs if exsit
    if 'additional_expr' in fitter[-1]:
        for k, v in list(fitter[-1]['additional_expr'].items()):
            params[k] = fit.Parameter(name=k, expr=v)

    ## do fitting
    r = _fitIQ(data, s, errors, params, func, silent=silent)
    r.add_functions(others)
    return r

def adjust_fitrange(nfwhm, ndata, nparams, peakparams, factor=1):
    """
    Adjust fit range to make sure :math:`N_{free} \geq nparam`

    :param nfwhm: a float, minimum fit length in unit of FWHM
    :param ndata: a integer, length of data
    :param nparams: a integer, number of parameter
    :param peakparams: a dict of peak information
    :param factor: a integer

    :Note: *specify 2* as **factor** if you fit complex variable (as 2d-vector)

    **peakparams** are:

    - **f0ind**: a integer, index of peak center
    - **dl**: a integer, index count to reach the left half-maximum
    - **dr**: a integer, index count to reach the right half-maximum

    :return: (rbegin, rend)

    - **rbegin** is a integer of fit range begin
    - **rend** is a integer of fit range end
    """
    f0ind = peakparams['f0ind']
    l, c, r = f0ind-peakparams['dl'], f0ind, f0ind+peakparams['dr']
    #l, c, r = f0ind-peakparams['dl']-1, f0ind, f0ind+peakparams['dr']+1
    if (r - l + 1)*factor >= nparams:
        rbegin, rend = int(c-nfwhm*(c-l)), int(c+nfwhm*(r-c))
    else:
        n = nparams/(float(r - l + 1)*factor)
        rbegin, rend = int(f0ind-nfwhm*n*(c-l)), int(c+nfwhm*n*(r-c))
    if rbegin < 0:
        if rend + (-rbegin) >= ndata:
            rend = ndata-1
            rbegin = 0
        else:
            rend = rend + (-rbegin)
            rbegin = 0
    if rend >= ndata:
        if rbegin - (rend - ndata) < 0:
           rbegin = 0
           rend = ndata-1
        else:
           rbegin = rbegin - (rend - ndata)
           rend = ndata-1

    if (rend - rbegin + 1)*factor < 11:
        print(f'ERROR:: {rbegin} to {rend} times {factor} < 11')
        raise MKIDDataException("Fit range guess error")

    return rbegin, rend

