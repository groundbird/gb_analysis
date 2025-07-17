import numpy as np
from .TODdata import TODdata
from .readfile import readfile_tod
from operator import attrgetter
from . import misc

def combine_TODdata(tod1,tod2):
    """
    combine two TODdata as the sequience (kind is fixed to `I-Q`)
    """
    if tod1 is None: return tod2
    if tod2 is None: return tod1
    if tod1.rate != tod2.rate or tod1.f != tod2.f:
        raise RuntimeError('Not same TOD data!')

    t1 = tod1.t
    i1 = tod1.i
    q1 = tod1.q

    t2 = tod2.t + tod1.t[-1] + 1./tod1.rate
    i2 = tod2.i
    q2 = tod2.q

    return TODdata(t1.tolist()+t2.tolist(),'I-Q',(i1.tolist()+i2.tolist(),q1.tolist()+q2.tolist()),tod1.f)

def _multreadfile_frq(kind, filenames, lo, **kws):

    if not hasattr(lo, '__iter__'):
        lo = [lo]*len(filenames)
    if len(lo) != len(filenames):
        raise RuntimeError(f'Different # of filenames ({len(filenames)}) and lo ({len(lo)}) in `_multreadfile`.')

    ret = list(zip(*[readfile_tod(kind,fname,-1,f,**kws) for f,fname in zip(lo,filenames)]))
    return [sorted(dlist,key=attrgetter('rate')) for dlist in ret] # sort with sampling rate

def _multreadfile_frqrate(kind, filenames, lo, **kws):

    ret = _multreadfile_frq(kind, filenames, lo, **kws)

    # duplication check with the first data
    ratelist = [d.rate for d in ret[0]]
    if len(ratelist) == len(set(ratelist)): return ret
    nrate = len(set(ratelist))

    cor_ret = [None] * len(ret)
    for i,dlist in enumerate(ret):
        cor_dlist = [None]*nrate
        for j,r in enumerate(set(ratelist)):
            rate_index = np.where([r==ir for ir in ratelist])
            combined = None
            for i in rate_index:
                combined = combine_TODdata(combined,dlist[i])
            cor_dlist[j] = combined
        cor_ret[i] = cor_dlist

    return cor_ret

def multreadfile_tod(kind, filenames, lo=None, merged=False, **kws):
    """
    read TODdata list from given files
    this function merges TOD with the same rate/KID-freq

    :param kind: one of ['rhea', 'riken', 'vna']
    :param filename: name of input file
    :param index: KID index number to be assigned in the input file
    """

    lo = misc.tofrq(lo)
    if not merged:
        return _multreadfile_frq(kind, filenames, lo, **kws)
    else:
        return _multreadfile_frqrate(kind, filenames, lo, **kws)

