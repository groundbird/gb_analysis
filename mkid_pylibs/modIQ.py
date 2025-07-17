import numpy as np
from . import misc

def subtract_blindtone2(mainfrq,mainIQ,lowfrq,lowIQ,highfrq,highIQ):
    blindIQ = lowIQ + (highIQ-lowIQ)/(highfrq-lowfrq)*(mainfrq-lowfrq)
    return subtract_blindtone(mainIQ,blindIQ)

def subtract_blindtone(mainIQ,blindIQ):
    #return _subtract_blindtone_JS(mainIQ,blindIQ)
    return _subtract_blindtone_HK(mainIQ,blindIQ)

def _subtract_blindtone_JS(mainIQ,blindIQ):
    pass

def _subtract_blindtone_HK(mainIQ,blindIQ):
    nz_on = 1.0/np.average(np.abs(mainIQ))*mainIQ
    nz_off = 1.0/np.average(np.abs(blindIQ))*blindIQ
    #del_anz_off = np.abs(nz_off) - np.average(np.abs(nz_off)[:100])
    del_anz_off = np.abs(nz_off) - np.average(np.abs(nz_off))
    anz_on_dash = np.abs(nz_on) - del_anz_off
    anz_on_dash2 = np.abs(nz_on) + del_anz_off

    #ctz_off = np.pi*np.sign(np.angle(blindIQ)) - np.angle(blindIQ)
    ctz_off = np.angle(blindIQ)
    del_ctz_off = ctz_off - np.average(ctz_off)
    #ctz_on = np.pi*np.sign(np.angle(mainIQ)) - np.angle(mainIQ)
    ctz_on = np.angle(mainIQ)
    ctz_on_dash = ctz_on - del_ctz_off
    ctz_on_dash2 = ctz_on + del_ctz_off

#    if np.std(anz_on_dash2) < np.std(anz_on_dash) and np.std(ctz_on_dash2) < np.std(ctz_on_dash):
#        anz_on_dash = anz_on_dash2
#        ctz_on_dash = ctz_on_dash2

    #tz_on_dash = -1*ctz_on_dash + np.array([0 if t>-np.pi else 2*np.pi for t in ctz_on_dash])
    az_on_dash = anz_on_dash*np.average(np.abs(mainIQ))
    tz_on_dash = ctz_on_dash

    # ctz_off = np.angle(blindIQ)
    # del_ctz_off = ctz_off - np.average(ctz_off[:100])
    # ctz_on = np.angle(mainIQ)
    # ctz_on_dash = ctz_on - del_ctz_off
    # tz_on_dash = ctz_on_dash
    # plt.plot(tz_on_dash,'skyblue')

    return az_on_dash*np.exp(1j*tz_on_dash)

def rwIQ(f,iq,swpft,deglitch=False):
    """
    add modified IQdata
    """
    return swpft.rewind(f, iq)

def modIQ(iq,ciq=None,f=None,cf=None):
    """
    add modified IQdata
    """
    if hasattr(cf, '__iter__'):
        lf = cf[0]
        liq = ciq[0]
        hf = cf[1]
        hiq = ciq[1]
        miq = subtract_blindtone2(f,iq,lf,liq,hf,hiq)
    else:
        miq = subtract_blindtone(iq,ciq)
    return miq

def add_moddata(f, data, swpft=None, cf=None, cdata=None, deglitch=False):
    """
    add modified IQdata
    - rwdata = rewound IQdata with swpft
    - mdata = 2ch subtracted IQdata with off-resonance
    - rwmdata = rewound subtracted IQdata

    :param f: list of trequency
    :param data: KidGlobalResponse to be modified
    :param swpft: KidFitResult for its swpIQdata
    :param cf: off-resonance frequency (or its list)
    :param cdata: off-resonance KidGlobalResponse (or its list)
    :param deglitch: The data will be deglitched bofre rewinding if `True`
    """

    is_fitted = (swpft is not None)
    is_cont   = (cdata is not None)

    if is_fitted:
        iq = data.iq
        if deglitch:
            amp,phs = misc.deglitch(data.x,[data.amplitude,data.phase])
            iq = amp*np.exp(1j*phs)
        rwiq = swpft.rewind(f, iq)
        # additional modification by H.Kutsuma
        #if 'circle_fit' in swpft.info:
        #rfit = swpft.info['circle_fit']['rfit']
        #zcfit = swpft.info['circle_fit']['zcfit']
        #rwiq = 1.0/rfit*(rwiq - zcfit)*np.exp(-1j*np.angle(zcfit))
        #rwiq = 1.0/rfit*(rwiq - rfit)*np.exp(-1j*np.angle(zcfit)) # fix translation to the origin
        data.add_IQdata(rwiq,'rwdata')

    if is_cont:
        if hasattr(cdata, '__iter__'):
            lf = cf[0]
            liq = cdata[0].iq
            hf = cf[1]
            hiq = cdata[1].iq
            miq = subtract_blindtone2(f,data.iq,lf,liq,hf,hiq)
        else:
            miq = subtract_blindtone(data.iq,cdata.iq)
        data.add_IQdata(miq,'mdata')

    if data.mdata is not None and is_fitted:
        iq = data.mdata.iq
        if deglitch:
            amp,phs = misc.deglitch(data.x,[data.mdata.amplitude,data.mdata.phase])
            iq = amp*np.exp(1j*phs)
        rwmiq = swpft.rewind(f,iq)
        data.add_IQdata(rwmiq,'rwmdata')
    pass

