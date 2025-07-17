# -*- coding: utf-8 -*-
#### fitter definitions
"""
Fitter fuction object

fitter named "ft" provides:
- ft_paramnames : a list of parameter names
- ft_paramlatex : a list of parameter names in latex notation
- ft_ewa = make_EWAft() : fit function expression
- ft_ewa_bg = make_EWAft_bg() : bkgd function expression
- ft_guess(data) : return a dictionary object with 'param_name:value'
- ft_rewind(x,y,...) : return rewound y
- fitter_ft : above functions are arranged in a `Fitter` namedtuple

available fitters: 'mazinrev', 'gao', 'gaolinbg', 'blank'
"""

import collections
import sympy
import numpy as np
from .fit.expr_with_args import Expr_with_args
from .peak_search import search_peak_center
from .peak_search import analytical_search_peak, analytical_search_peak2
from .misc import circle_fit

Fitter = collections.namedtuple('Fitter', 'func, guess, paramnames, info')

complex_fitters = ['gaolinbg2f', 'gaolinbg2l', 'gaolinbg', 'blank', 'mazinrev', 'gao']
real_fitters    = []
all_fitters     = complex_fitters + real_fitters

#################################################################
# gaolinbg2(f or l): a complex function from Gao's D thesis, plus linear term for 2 KID made by ysueno.
################################################################
gaolinbg2_paramnames = 'arga absa tau fr1 Qr1 Qc1 phi01 fr2 Qr2 Qc2 phi02 c'.split()
gaolinbg2_paramlatex = r'(\arg{a}) |a| tau f_r1 Q_r1 Q_c1 phi_01 f_r2 Q_r2 Q_c2 phi_02 c'

def make_EWA_gaolinbg2f():
    arg_names = gaolinbg2_paramnames
    arg_latex = gaolinbg2_paramlatex

    ## define symbols to use as parameters
    arg_symbols = sympy.symbols(arg_latex)
    arga, absa, tau, fr1, Qr1, Qc1, phi01, fr2, Qr2, Qc2, phi02, c = arg_symbols

    from sympy import exp, I, pi

    x = sympy.symbols('x')
    expr = (absa * exp(-I*(2*pi*x*tau - arga))* 
            (1+c*(x-fr1)-Qr1/Qc1*exp(I*phi01)/(1+2*I*Qr1*((x-fr1)/fr1))
            - Qr2/Qc2*exp(I*phi02)/(1+2*I*Qr2*((x-fr2)/fr2)))
    )

    return Expr_with_args(expr, arg_symbols, arg_names)

def make_EWA_gaolinbg2l():
    arg_names = gaolinbg2_paramnames
    arg_latex = gaolinbg2_paramlatex

    ## define symbols to use as parameters
    arg_symbols = sympy.symbols(arg_latex)
    arga, absa, tau, fr1, Qr1, Qc1, phi01, fr2, Qr2, Qc2, phi02, c = arg_symbols

    from sympy import exp, I, pi

    x = sympy.symbols('x')
    expr = (absa * exp(-I*(2*pi*x*tau - arga))* 
            (1+c*(x-fr2)-Qr1/Qc1*exp(I*phi01)/(1+2*I*Qr1*((x-fr1)/fr1))
            - Qr2/Qc2*exp(I*phi02)/(1+2*I*Qr2*((x-fr2)/fr2)))
    )

    return Expr_with_args(expr, arg_symbols, arg_names)

def make_EWABG_gaolinbg2():
    arg_names = gaolinbg2_paramnames
    arg_latex = gaolinbg2_paramlatex

    ## define symbols to use as parameters
    arg_symbols = sympy.symbols(arg_latex)
    arga, absa, tau, fr1, Qr1, Qc1, phi01, fr2, Qr2, Qc2, phi02, c = arg_symbols

    from sympy import exp, I, pi

    x = sympy.symbols('x')
    expr = (absa * exp(-I*(2*pi*x*tau - arga))*(1+c*(x-fr1)))

    return Expr_with_args(expr, arg_symbols, arg_names)

def guess_gaolinbg2f(data, dep = 3):
    print(f'dep : {dep}')
    x       = data.x
    deltax  = x[1] - x[0]
    #pdict   = search_peak_center(data)
    pdict   = analytical_search_peak2(data, dep = dep)
    y1      = data.iq[pdict['f1ind']]
    FWHM1    = pdict['f1']/pdict['Q1']
    FWHM2    = pdict['f2']/pdict['Q2']
    ddeg    = data.deg[1:] - data.deg[:-1]
    tau     = -np.average(ddeg[abs(ddeg)<180])*np.pi/180.0/deltax/2/np.pi
    f1      = pdict['f1']
    f2      = pdict['f2']
    theta   = np.angle(y1)
    arga    = np.angle(y1*np.exp(1j*2*np.pi*tau*f1))
    absa    = pdict['a_off']
    fr1     = f1
    fr2     = f2
    Qr1     = f1/FWHM1
    Qr2     = f2/FWHM2
    Qc1     = Qr1
    Qc2     = Qr2
    phi01   = 0
    phi02   = 0
    c       = 0

    return dict(list(zip(gaolinbg2_paramnames, (arga, absa, tau, fr1, Qr1, Qc1, phi01, fr2, Qr2, Qc2, phi02, c))))

def guess_gaolinbg2l(data, dep = 3):
    x       = data.x
    deltax  = x[1] - x[0]
    #pdict   = search_peak_center(data)
    pdict   = analytical_search_peak2(data, dep = dep)
    y2      = data.iq[pdict['f2ind']]
    FWHM1    = pdict['f1']/pdict['Q1']
    FWHM2    = pdict['f2']/pdict['Q2']
    ddeg    = data.deg[1:] - data.deg[:-1]
    tau     = -np.average(ddeg[abs(ddeg)<180])*np.pi/180.0/deltax/2/np.pi
    f1      = pdict['f1']
    f2      = pdict['f2']
    theta   = np.angle(y2)
    arga    = np.angle(y2*np.exp(1j*2*np.pi*tau*f2))
    absa    = pdict['a_off']
    fr1     = f1
    fr2     = f2
    Qr1     = f1/FWHM1
    Qr2     = f2/FWHM2
    Qc1     = Qr1
    Qc2     = Qr2
    phi01   = 0
    phi02   = 0
    c       = 0

    return dict(list(zip(gaolinbg2_paramnames, (arga, absa, tau, fr1, Qr1, Qc1, phi01, fr2, Qr2, Qc2, phi02, c))))

def rewind_gaolinbg2f(x, y, arga, absa, tau, fr1, Qr1, Qc1, phi01, fr2, Qr2, Qc2, phi02, c):
    tmp = y/absa/np.exp(-1j*(2*np.pi*x*tau-arga)) - c*(x-fr1) + Qr2/Qc2*np.exp(1j*phi02)/(1+2*1j*Qr2*((x-fr2)/fr2))
    return (tmp-1)*Qc1/Qr1/np.exp(1j*phi01) + 0.5

def rewind_gaolinbg2l(x, y, arga, absa, tau, fr1, Qr1, Qc1, phi01, fr2, Qr2, Qc2, phi02, c):
    tmp = y/absa/np.exp(-1j*(2*np.pi*x*tau-arga)) - c*(x-fr2) + Qr1/Qc1*np.exp(1j*phi01)/(1+2*1j*Qr1*((x-fr1)/fr1))
    return (tmp-1)*Qc2/Qr2/np.exp(1j*phi02) + 0.5

ewa_gaolinbg2f    = make_EWA_gaolinbg2f()
ewa_gaolinbg2l    = make_EWA_gaolinbg2l()
ewabg_gaolinbg2  = make_EWABG_gaolinbg2()

# focus on first KID
fitter_gaolinbg2f = Fitter(ewa_gaolinbg2f, guess_gaolinbg2f, gaolinbg2_paramnames,
                         {'bgfunc'         : ewabg_gaolinbg2,
                          'rewindfunc'     : rewind_gaolinbg2f,
                          'additional_expr': {'Qi1': '1/(1/Qr1 - 1/Qc1*cos(phi01))'},
                          'prefit_and_fix': []
                          })

# focus on last KID
fitter_gaolinbg2l = Fitter(ewa_gaolinbg2l, guess_gaolinbg2l, gaolinbg2_paramnames,
                         {'bgfunc'         : ewabg_gaolinbg2,
                          'rewindfunc'     : rewind_gaolinbg2l,
                          'additional_expr': {'Qi2': '1/(1/Qr2 - 1/Qc2*cos(phi02))'},
                          'prefit_and_fix': []
                          })

#################################################################
# gaolinbg: a complex function from Gao's D thesis, plus linear term
################################################################
gaolinbg_paramnames = 'arga absa tau fr Qr Qc phi0 c'.split()
gaolinbg_paramlatex = r'(\arg{a}) |a| tau f_r Q_r Q_c phi_0 c'

def make_EWA_gaolinbg():
    arg_names = gaolinbg_paramnames
    arg_latex = gaolinbg_paramlatex

    ## define symbols to use as parameters
    arg_symbols = sympy.symbols(arg_latex)
    arga, absa, tau, fr, Qr, Qc, phi0, c = arg_symbols

    from sympy import exp, I, pi

    x = sympy.symbols('x')
    expr = (absa * exp(-I*(2*pi*x*tau - arga))*
            (1+c*(x-fr)-Qr/Qc*exp(I*phi0)/(1+2*I*Qr*((x-fr)/fr))))

    return Expr_with_args(expr, arg_symbols, arg_names)

def make_EWABG_gaolinbg():
    arg_names = gaolinbg_paramnames
    arg_latex = gaolinbg_paramlatex

    ## define symbols to use as parameters
    arg_symbols = sympy.symbols(arg_latex)
    arga, absa, tau, fr, Qr, Qc, phi0, c = arg_symbols

    from sympy import exp, I, pi

    x = sympy.symbols('x')
    expr = (absa * exp(-I*(2*pi*x*tau - arga))*(1+c*(x-fr)))

    return Expr_with_args(expr, arg_symbols, arg_names)

def guess_gaolinbg(data):
    x       = data.x
    deltax  = x[1] - x[0]
    #pdict   = search_peak_center(data)
    pdict   = analytical_search_peak(data)
    y0      = data.iq[pdict['f0ind']]
    FWHM    = pdict['f0']/pdict['Q']
    ddeg    = data.deg[1:] - data.deg[:-1]
    tau     = -np.average(ddeg[abs(ddeg)<180])*np.pi/180.0/deltax/2/np.pi
    f0      = pdict['f0']
    theta   = np.angle(y0)
    arga    = np.angle(y0*np.exp(1j*2*np.pi*tau*f0))
    absa    = pdict['a_off']
    fr      = f0
    Qr      = f0/FWHM
    Qc      = Qr
    phi0    = 0
    c       = 0

    return dict(list(zip(gaolinbg_paramnames, (arga, absa, tau, fr, Qr, Qc, phi0, c))))

def rewind_gaolinbg(x, y, arga, absa, tau, fr, Qr, Qc, phi0, c):
    tmp = y/absa/np.exp(-1j*(2*np.pi*x*tau-arga)) - c*(x-fr)
    return (tmp-1)*Qc/Qr/np.exp(1j*phi0) + 0.5

ewa_gaolinbg    = make_EWA_gaolinbg()
ewabg_gaolinbg  = make_EWABG_gaolinbg()

fitter_gaolinbg = Fitter(ewa_gaolinbg, guess_gaolinbg, gaolinbg_paramnames,
                         {'bgfunc'         : ewabg_gaolinbg,
                          'rewindfunc'     : rewind_gaolinbg,
                          'additional_expr': {'Qi': '1/(1/Qr - 1/Qc*cos(phi0))'},
                          'prefit_and_fix': []
                          })

################################################################
# blank: a complex function to fit baseline for cable delay
################################################################
blank_paramnames = 'arga absa tau c'.split()
blank_paramlatex = r'(\arg{a}) |a| tau c'
def make_EWA_blank():
    arg_names = blank_paramnames
    arg_latex = blank_paramlatex

    ## define symbols to use as parameters
    arg_symbols = sympy.symbols(arg_latex)
    arga, absa, tau, c = arg_symbols

    from sympy import exp, I, pi

    x = sympy.symbols('x')
    expr = (absa * exp(-I*(2*pi*x*tau - arga))*(1+c*x))

    return Expr_with_args(expr, arg_symbols, arg_names)

def guess_blank(data):
    x       = data.x
    deltax  = x[1] - x[0]

    y = data.amplitude
    a0 = (y[-1] - y[0])/(x[-1] - x[0])
    b0 = y[0] - a0*x[0]

    f0ind = np.argmin(y - (a0*x+b0))
    f0 = x[f0ind]
    y0 = y[f0ind]
    ddeg    = data.deg[1:] - data.deg[:-1]
    tau     = -np.average(ddeg[abs(ddeg)<180])*np.pi/180.0/deltax/2/np.pi
    arga    = np.angle(y0*np.exp(1j*2*np.pi*tau*f0))
    absa    = b0
    c       = a0/absa

    return dict(list(zip(blank_paramnames, (arga, absa, tau, c))))

def rewind_blank(x, y, arga, absa, tau, c):
    tmp = y/absa/np.exp(-1j*(2*np.pi*x*tau-arga))/(1+c*x)
    return tmp

ewa_blank    = make_EWA_blank()

fitter_blank = Fitter(ewa_blank, guess_blank, blank_paramnames,
                      {'rewindfunc': rewind_blank,
                       # 'bgfunc'         : ,
                       # 'additional_expr': {},
                       # 'prefit_and_fix': []
                       })


#################################################################
#################################################################
############ NOT SUPPORTED BELOW FUNCTIONS SO FAR ###############
#################################################################
#################################################################


#################################################################
# mazinrev: a complex function from Mazin's D thesis
################################################################
mazinrev_paramnames = ['FWHM', 'f0', 'a_on', 'a_off', 'v', 'c', 'theta', 'gr', 'Ic', 'Qc']
mazinrev_paramlatex = r'FWHM, f_0 a_{on} a_{off} v c theta g_r I_c Q_c'
def make_EWAmazinrev():
    arg_names = mazinrev_paramnames
    arg_latex = mazinrev_paramlatex

    ## define symbols to use as parameters
    arg_symbols = sympy.symbols(arg_latex)
    FWHM, f0, a_on, a_off, v, c, theta, gr, Ic, Qc = arg_symbols

    from sympy import exp, I, pi, re, im

    x = sympy.symbols('x')
    deltax = x - f0
    origx = f0
    w = deltax/FWHM
    f = a_off + (a_off - a_on)*(2*I*w/(1+2*I*w) - 1) + c*deltax
    expr = (re(f)+I*gr*im(f))*exp(I*(theta-v*deltax)) + (Ic + I*Qc)
    return Expr_with_args(expr, arg_symbols, arg_names)

def make_EWAmazinrev_bg():
    arg_latex = mazinrev_paramlatex
    arg_names = mazinrev_paramnames

    ## define symbols to use as parameters
    arg_symbols = sympy.symbols(arg_latex)
    FWHM, f0, a_on, a_off, v, c, theta, gr, Ic, Qc = arg_symbols

    from sympy import exp, I, pi, re, im

    x = sympy.symbols('x')
    deltax = x - f0
    origx = f0
    w = deltax/FWHM
    f = a_off + c*deltax
    expr = (re(f)+I*gr*im(f))*exp(I*(theta-v*deltax)) + (Ic + I*Qc)
    return Expr_with_args(expr, arg_symbols, arg_names)

def mazinrev_guess(data):
    x      = data.x
    deltax = x[1] - x[0]
    pdict   = search_peak_center(data)
    y0      = data.iq[pdict['f0ind']]
    FWHM    = pdict['f0']/pdict['Q']
    ddeg    = data.deg[1:] - data.deg[:-1]
    v       = -np.average(ddeg[abs(ddeg)<180])*np.pi/180.0/deltax
    theta   = np.angle(y0)
    gr      = 1.0
    y_unrev = data.iq * np.exp(-1j*(-v*(x-pdict['f0']) + theta))
    c       = 0.0
    Ic, Qc  = 0.0, 0.0
    return dict(list(zip(mazinrev_paramnames, (FWHM, pdict['f0'], pdict['a_on'], pdict['a_off'], v, c, theta, gr, Ic, Qc))))

def mazinrev_rewind(x, y, FWHM, f0, a_on, a_off, v, c, theta, gr, Ic, Qc):
    deltax = x - f0
    f2  = (y - (Ic + 1j*Qc))/np.exp(1j*(theta-v*deltax))
    f   = np.real(f2) + 1j*np.imag(f2)/gr
    y_  = f - c*deltax
    return (y_ - a_off)/(a_off - a_on) + 0.5

mazinrev_ewa    = make_EWAmazinrev()
mazinrev_ewa_bg = make_EWAmazinrev_bg()
fitter_mazinrev = Fitter(mazinrev_ewa, mazinrev_guess, mazinrev_paramnames,
                         {'bgfunc': mazinrev_ewa_bg, 'rewindfunc': mazinrev_rewind})

#################################################################
# gao: a complex function from Gao's D thesis
################################################################
gao_paramnames = 'arga absa tau fr Qr Qc phi0'.split()
gao_paramlatex = r'(\arg{a}) |a| tau f_r Q_r Q_c phi_0'

def make_EWAgao():
    arg_names = gao_paramnames
    arg_latex = gao_paramlatex

    ## define symbols to use as parameters
    arg_symbols = sympy.symbols(arg_latex)
    arga, absa, tau, fr, Qr, Qc, phi0 = arg_symbols

    from sympy import exp, I, pi

    x = sympy.symbols('x')
    expr = (absa * exp(-I*(2*pi*x*tau - arga))*
            (1-Qr/Qc*exp(I*phi0)/(1+2*I*Qr*((x-fr)/fr))))

    return Expr_with_args(expr, arg_symbols, arg_names)

def make_EWAgao_bg():
    arg_names = gao_paramnames
    arg_latex = gao_paramlatex

    ## define symbols to use as parameters
    arg_symbols = sympy.symbols(arg_latex)
    arga, absa, tau, fr, Qr, Qc, phi0 = arg_symbols

    from sympy import exp, I, pi

    x = sympy.symbols('x')
    expr = (absa * exp(-I*(2*pi*x*tau - arga)))

    return Expr_with_args(expr, arg_symbols, arg_names)

def gao_guess(data):
    x      = data.x
    deltax = x[1] - x[0]
    pdict   = search_peak_center(data)
    y0      = data.iq[pdict['f0ind']]
    FWHM    = pdict['f0']/pdict['Q']
    ddeg    = data.deg[1:] - data.deg[:-1]
    tau     = -np.average(ddeg[abs(ddeg)<180])*np.pi/180.0/deltax/2/np.pi
    f0      = pdict['f0']
    theta   = np.angle(y0)
    arga    = np.angle(y0*np.exp(1j*2*np.pi*tau*f0))
    absa    = pdict['a_off']
    fr      = f0
    Qr      = f0/FWHM
    Qc      = Qr
    phi0    = 0

    return dict(list(zip(gao_paramnames, (arga, absa, tau, fr, Qr, Qc, phi0))))

def gao_rewind(x, y, arga, absa, tau, fr, Qr, Qc, phi0):
    tmp = y/absa/np.exp(-1j*(2*np.pi*x*tau-arga))
    return (tmp-1)*Qc/Qr + 0.5

gao_ewa    = make_EWAgao()
gao_bgewa  = make_EWAgao_bg()
fitter_gao = Fitter(gao_ewa, gao_guess, gao_paramnames,
                    {'bgfunc': gao_bgewa, 'rewindfunc': gao_rewind,
                     'additional_expr': {'Qi': '1/(1/Qr - 1/Qc*cos(phi0))'},
                     'prefit_and_fix': [],
                    })

