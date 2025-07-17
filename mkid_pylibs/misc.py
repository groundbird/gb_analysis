# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import exceptions
import builtins as exceptions
from scipy.constants import physical_constants
from scipy.optimize import leastsq

#from . import fit

#### python tool

def mult_getattr(x,d):
    if d == '': return x
    tmpary = d.split('.')
    cur_attr = tmpary[0]
    next_attr = '.'.join(tmpary[1:]) if len(tmpary)>0 else ''
    return mult_getattr(getattr(x,cur_attr),next_attr)

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

def down_sample(x, nsample):
    if hasattr(x, 'down_sample'):
        return x.down_sample(nsample)
    else:
        D = []
        for i in range(int(np.ceil(len(x)/float(nsample)))):
            beg = i*nsample
            end = min(len(x)-1, (i+1)*nsample)
            D.append(np.average(x[beg:end]))
        return np.array(D)

def deglitch(xs, yss, sources=None,
             baseline_thresh=6.0, glitch_thresh=5.0, clusterize_thresh=2):
    """
    Deglitch `yss` using `sources`, assuming there are glitches
    at the same time of all data in `yss`.
    Too close glitches or broad glitch are treated as one glitch.

    :param yss: an array of 1-D arrays of data to deglitch
    :param sources: source data to detect glitch. If None, yss is used
    :param baseline_thresh: threshold to decide baseline (`baseline_threshold`*sigma)
    :param glitch_thresh: threshold to decide as glitch (`glitch_threshold`*sigma)
    :param clusterize_thresh: if gap between glitches are less than or equal to this,
                              treat them as one glitch.

    Return:
        results: an array of 1-D arrays, that is deglitched yss
    """

    if sources is None:
        sources = yss

    def _numdif_2(y, dx):
        return (y[2:]+y[:-2]-2*y[1:-1])/dx**2

    ave = np.average(sources, axis=0)
    dx  = xs[1] - xs[0]
    diff2 = np.array(_numdif_2(ave, dx))
    sigma = np.std(diff2)
    good  = np.abs(diff2) < (baseline_thresh*sigma)
    sigma = np.std(diff2[good])
    bad   = (np.abs(diff2) >= glitch_thresh*sigma)
    ## treat broad glitch (or too close glitches) as one glitch
    bad   = misc.clusterizer(bad, clusterize_thresh)
    good  = np.logical_not(bad)
    results = []

    for ys in yss:
        center = np.interp(xs[1:-1], xs[1:-1][good], ys[1:-1][good])
        deglitched = np.concatenate(([ys[0]], center, [ys[-1]]))
        results.append(deglitched)

    return results

def oddify(n):
    return int(np.ceil(n/2.0)*2+1)

def ftoi(ff):
    return int(np.floor(ff+0.0001))

### Booleans

def widener(bool_array):
    '''
    Convert bool_array from
    [False, ... , False, False,  True, ... ,  True, False, False, ... , False]
    to
    [False, ... , False,  True,  True, ... ,  True,  True, False, ... , False]
    '''
    new_array = bool_array.copy()
    for i in range(len(bool_array)-1):
        if (bool_array[i] == False) & (bool_array[i+1] == True):
            new_array[i] = True
        elif (bool_array[i] == True) & (bool_array[i+1] == False):
            new_array[i+1] = True
        else:
            pass
    return new_array

def clusterizer(bool_array,threshold=1):
    """
    Fill `False` small gap (within `threshold`) with `True` in bol_array
    e.g.
    [True, False, False, True, True]
    to
    [True, False, False, True, True] (if threshold == 1)
    [True, True, True, True, True] (if threshold >= 2)

    :param threshold: an interger, allowed gap between True's in array
    """

    new_array = np.copy(bool_array)
    prev  = 0
    first = True
    for i, x in enumerate(new_array):
        if x:
            if (not first) and i - prev <= threshold + 1:
                for j in range(prev, i):
                    new_array[j] = True
            prev  = i
            first = False
    return new_array

#### KID analysis tool

def linearized_phase(x):
    return 2*np.tan(x/2.)

def amplitude_to_dB(x):
    return 20.0*np.log(x)/np.log(10.0)

def dB_to_amplitude(x):
    return 10.0**(x/20.0)

def get_data_type(IQdata,kind=None):
    '''
    Return specified IQdata.
    :param IQdata: The Swpdata/TODdata to be treated.
    :param kind: one of ['rwmdata', 'rwdata', 'mdata', 'data']
                 'rw': rewound one
                 'm' : modified by using off-resonance
                 if not specified, the priority is 'rwmdata' --> 'rwdata' --> 'mdata' --> 'data'
    '''
    data = None
    if (kind is None and data is None) or kind == 'rwmdata':
        data = IQdata.rwmdata
        if data is not None: kind='rwmdata'
        pass
    if (kind is None and data is None) or kind == 'rwdata':
        data = IQdata.rwdata
        if data is not None: kind='rwdata'
        pass
    if (kind is None and data is None) or kind == 'mdata':
        data = IQdata.mdata
        if data is not None: kind='mdata'
        pass
    if (kind is None and data is None) or kind == 'data' or kind == 'rawdata':
        data = IQdata
        if data is not None: kind='data'
        pass
    if data is None:
        try:
            if kind == 'fitdata.rwdata':
                data = IQdata.fitdata.rwdata
                if data is not None: kind='fitdata.rwdata'
                pass
            if kind == 'fitdata' or kind == 'fitdata.data' or kind == 'fitdata.rawdata':
                data = IQdata.fitdata
                if data is not None: kind='fitdata'
                pass
        except Exception as e:
            print(e)
            pass
    if data is None:
        raise RuntimeError(f'ERROR: Invalid/Not-assigned IQdata type: kind = {kind} for {IQdata.__class__.__name__}')

    return data,kind

def get_datatype_description(kind):
    if kind == 'data':
        return 'raw'
    ss = ''
    if 'rw' in kind:
        ss += 'rewind,'
    if 'm' in kind:
        ss += 'mod,'
    if 'fit' in kind:
        ss += 'fit,'
    return ss[:-1]

def convert_Kelvin_nqp(T,material='Al'): # T[K]
    kB = physical_constants['Boltzmann constant in eV/K'][0]
    if material == 'Al':
        N0 = 1.74e10 #[/eV/um^3]
        D0 = 193e-6 #[eV]
    else:
        raise('ERROR:: no materials other than Al are not supported yet. add your material ({material}) in misc.py.')

    return 2*N0*np.sqrt(2*np.pi*kB*T*D0)*np.exp(-D0/(kB*T)) #[/um^3]

def convert_tqp_nqp(tqp,material='Al'): # tqp[sec]
    kB = physical_constants['Boltzmann constant in eV/K'][0]
    if material == 'Al':
        N0 = 1.74e10 #[/eV/um^3]
        D0 = 193e-6 #[eV]
        Tc = 1.2 #[K]
        t0 = 450e-9 #[sec]
    else:
        raise('ERROR:: no materials other than Al are not supported yet. add your material ({material}) in misc.py.')

    return (t0/tqp)*N0*((kB*Tc)**3)/(2*(D0**2)) #[/um^3]

def convert_Nqp_P(Nqp,tqp,material='Al'): # Nqp[], tqp[sec]
    if material == 'Al':
        D0 = 2.898e-23 #[J]
        eta = 0.4 #[]
    else:
        raise('ERROR:: no materials other than Al are not supported yet. add your material ({material}) in misc.py.')

    return Nqp*D0/(eta*tqp) #[J/s]=[W]

def circle_fit(data):
    fitfunc = lambda p, x, y: p[0]*x + p[1]*y + p[2]
    errfunc = lambda p, x, y: fitfunc(p, x, y) - (x**2 + y**2)
    A_ini = np.max(data.real) + np.min(data.real)
    B_ini = np.max(data.imag) + np.min(data.imag)
    C_ini = ((np.max(data.real)-np.min(data.real))/2.0)**2 - A_ini**2 - B_ini**2
    p0 = [A_ini, B_ini, C_ini]
    A, B, C = leastsq(errfunc, p0[:], args=(data.real, data.imag))[0]
    zc = A/2.0 + 1j* B/2.0
    r = np.sqrt(C + (A/2.0)**2 + (B/2.0)**2)
    return zc, r

def linear_fit(xs, ys, exs, eys, xname='x', yname='y',
               silent=False, fitname='lin', plotted=False):
    '''
    Fit linear function.
    :param xs: x variable to be fitted
    :param ys: y variable to be fitted
    :param exs: x error to be fitted
    :param eys: y error to be fitted
    :param xname: x variable name (used if plotted is True)
    :param yname: y variable name (used if plotted is True)
    :param silent: If False, detail fitting information will be dumped.
    :param fitname: One of [lin, c+lin]
    :param plotted: If not False, fitting result will be plotted, and saved as [plotted].pdf.
    '''

    from . import fit
    ## Initial values and fitting function
    i_aa = (ys[0]-ys[-1])/(xs[0]-xs[-1])
    i_bb = ys[0] - i_aa*xs[0]
    fitfunc = None
    if fitname == 'lin':
        def fitfunc(xs, a=i_aa, b=i_bb):
            fys = a*xs+b
            return fys
    elif fitname == 'c+lin':
        i_cc = np.mean(xs)
        def fitfunc(xs, a=i_aa, b=i_bb, c=i_cc):
            fys = (a*xs+b) * (xs>c) + (a*c+b) * (xs<=c)
            return fys
    else:
        raise RuntimeError(f'Invalid fitname: {fitname}. fitname option should be one of [lin, c+lin].')

    ## Fitting
    result = fit.fit(fitfunc,xs,ys,err=eys,silent=silent)

    ## Plot
    if plotted:
        fig, axs = plt.subplots()
        fig.set_size_inches(10,6)
        axs.errorbar(xs,ys,xerr=exs,yerr=eys,fmt='bo',ecolor='b')
        fitx = np.arange(xs.min(),xs.max(),np.abs(xs.min()-xs.max())/400)
        axs.plot(fitx,result.eval(fitx),color='r',marker='')
        axs.set_xlabel(xname)
        axs.set_ylabel(yname)
        axs.grid(ls=':')
        fig.tight_layout()
        return result, fig, axs

    return result

def is_difflarge(ys1,ys2,nsigma=5):
    diff = np.array(ys1)-np.array(ys2)
    return np.abs(diff - diff.mean()) < nsigma*diff.std()

#### frequency tool

_unitdict = dict()
_unitdict['G'] = 1e9
_unitdict['M'] = 1e6
_unitdict['k'] = 1e3
_unitdict['m'] = 1e-3
_unitdict['u'] = 1e-6
_unitdict['n'] = 1e-9
def convunit(ss):
    if ss not in _unitdict:
        import warnings
        warnings.warn(f'invalid unit: {ss}. Choose from {_unitdict.keys()}')
        return 1
    return _unitdict[ss]
def tofrq(ss):
    if type(ss) is str:
        if ss[-2:].lower() == 'hz':
            ss = ss[:-2]
        if ss[-1] in _unitdict:
            return float(ss[:-1])*_unitdict[ss[-1]]
    return float(ss)

def tosi(d, sep='', ret_int = True):
    '''
    Convert number to string with SI prefix
    '''
    incPrefixes = ['k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']
    decPrefixes = ['m', 'µ', 'n', 'p', 'f', 'a', 'z', 'y']

    if d==0: return str(0)

    degree = int(np.floor(np.log10(np.abs(d)) / 3))

    scaled = d
    prefix = ''

    if degree != 0:
        ds = degree/np.fabs(degree)
        if ds == 1:
            if degree - 1 < len(incPrefixes):
                prefix = incPrefixes[degree - 1]
            else:
                prefix = incPrefixes[-1]
                degree = len(incPrefixes)
        elif ds == -1:
            if -degree - 1 < len(decPrefixes):
                prefix = decPrefixes[-degree - 1]
            else:
                prefix = decPrefixes[-1]
                degree = -len(decPrefixes)

        scaled = float(d * np.power(float(1000), -degree))
        prefix = sep + prefix
        pass

    if ret_int: scaled = int(scaled)
    s = f'{scaled}{sep}{prefix}'

    return(s)

#### debug tool

def print_line():
    print('==============================================================================')

#### figure tool

def def_style():
    import matplotlib as mpl
    mpl.rcParams.update({'font.size': 14})
    mpl.rcParams.update({'axes.facecolor':'w'})
    mpl.rcParams.update({'axes.edgecolor':'k'})
    mpl.rcParams.update({'figure.facecolor':'w'})
    mpl.rcParams.update({'figure.edgecolor':'w'})
    mpl.rcParams.update({'axes.grid':True})
    mpl.rcParams.update({'grid.linestyle':':'})
    mpl.rcParams.update({'figure.figsize':[12,9]})

def make_figure(fsz=(8,6),nc=1,nr=1,dpi=100):
    fig, axs = plt.subplots(facecolor='w', edgecolor='k', dpi=dpi, figsize=fsz, ncols=nc, nrows=nr)
    if nc>1 or nr>1:
        axs = axs.flatten()
    return fig, axs

def adj_axis(fig,xlabel=None,ylabe=None):
    pass

def adj_figure(fig,title='',tightlayout=True):#,alignedlabel=True):
    if fig is None:
        import warnings
        warnings.warn('figure is None')
        return
    if tightlayout:
        fig.tight_layout()
#    if alignedlabel:
#        fig.align_ylabels()
    if title !='':
        fig.suptitle(title,fontsize=16)
        hh = fig.get_figheight()
        fig.subplots_adjust(top=1-50/hh*0.01)

def show_figure_iterm(fig):
    try:
        fig.savefig(".tmp.png", bbox_inches='tight', pad_inches=0)
        from imgcat import imgcat
        imgcat(open(".tmp.png"))
    except Exception:
        import warnings
        warnings.warn('cannot call imgcat in your environment')
        pass

#### exception class
class MKIDDataException(exceptions.Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

#### test functions
def guess_cabledelay(sd):
    """
    Guess cable delay v from SweepData
    """
    dx     = sd.x[1] - sd.x[0]
    dangle = sd.deg[1:] - sd.deg[:-1]
    return -np.average(dangle[abs(dangle)<180])*np.pi/180.0/dx

def plot_iq(iq, *args, **kwargs):
    plt.plot(np.real(iq), np.imag(iq), *args, **kwargs)

def get_dtheta(r, fx):
    theta_r  = np.angle(r.fitted(r.fitparams['fr'].value))
    theta_fx = np.average(np.angle(fx.iq))

    return theta_r - theta_fx

def patch_if(data, cond):
    indices = np.argwhere(cond)
    result = data[:]
    for i in indices:
        if i == 0:
            result[i] = result[i+1]
        elif i == len(data)-1:
            result[i] = result[i-1]
        else:
            result[i] = (data[i+1]+data[i-1])/2.0
    return result

class ListTable(list):
    """
    Overridden list class which takes a 2-dimensional list of the form [[1,2,3],[4,5,6]], and renders an HTML Table in IPython Notebook.
    """

    def _repr_html_(self):
        html = ["<table>"]
        for row in self:
            html.append("<tr>")

            for col in row:
                html.append("<td>{0}</td>".format(col))

            html.append("</tr>")
        html.append("</table>")
        return ''.join(html)

def summary(rs):
    res = ListTable()
    res.append( ('fr', 'Qr', 'Qi', 'phi0', 'WSSR') )
    import ad.admath as am
    for r in rs:
        c = ErrorCalculator(r.fitparams)
        Qi = 1/(1/c['Qr'] - 1/c['Qc']*am.cos(c['phi0']))
        res.append((c.str(c['fr']),c.str(c['Qr']),c.str(Qi),c.str(c['phi0']), r.info['s_sq']))
    return res

def rebin_log(x,y, begin=None, base=2, factor=2, average=True):
    """
    Rebin x, y (x is isospaced)

    :param begin: at what x rebin starts
    :param base: rebin bin width changes at begin*(base**n), n=0, 1, 2, ...
    :param factor: rebin bin width multiplier
    :param average: do average over new bin?
    """
    x_, y_ = [], []
    ind = 0
    width = 1
    if begin is None:
        nextx = 10**(np.floor(np.log(x[0])/np.log(10.0))+2)
    else:
        nextx = begin
    while ind < len(y):
        if x[ind] >= nextx:
            width *= factor
            nextx = nextx*base
        bintop = min(ind+width, len(y))
        binwidth = bintop-ind
        x_.append(x[ind])
        if average:
            y_.append(sum(y[ind:bintop]/float(binwidth)))
        else:
            y_.append(sum(y[ind:bintop]))
        ind += binwidth
    return np.array(x_), np.array(y_)
