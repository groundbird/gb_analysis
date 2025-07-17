import numpy as np
from . import misc
from . import fit

def fit_trig(toddata, kind=None, i_tqp=100e-6, silent=False, linearized=False, **kws):
    '''
    THIS IS FOR TRIGGERED TOD!! NOT USE WITH USUAL TOD DATA.
    Fit one triggered phase TOD data.
    :param toddata: TodData to be treated.
    :param kind: TodData type to be calculated.
    :param i_tqp: initial value for the quasi-particle life time for fitting.
    :param silent: If False, detail fitting information will be dumped.
    :param linearized: If True, the phase will be linearized with misc.linearized_phase().
    '''

    trigdata,kind = misc.get_data_type(toddata,kind=kind)
    toddata.info['type_trigfit'] = kind

    ## data to be fitted
    xs = toddata.time #[sec]
    if linearized:
        ys = trigdata.lincorphase
    else:
        ys = trigdata.corphase

    ## Initial values and fitting function
    def fitfunc(x, t0=100*1e-6, tqp=i_tqp, n0=np.max(ys)-(ys).mean(), bkgd=(ys).mean()):
        return np.heaviside(x-t0,1)*(n0*np.exp(-(x-t0)/tqp)) + bkgd

    r = fit.fit(fitfunc, xs, ys, silent=silent, **kws)

    return r

