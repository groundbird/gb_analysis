import numpy as np
from . import misc
from .PSDdata import PSDdata

def calc_nep(psd, dfr_dNqp, dQiInv_dNqp, swpfitresult, volume=None, T=None, tqp=None, Nqp=None):

    Qi = swpfitresult.params['Qi'].value
    fr = swpfitresult.params['fr'].value
    tres = Qi/(np.pi*fr) #[sec]

    P = misc.convert_Nqp_P(Nqp,tqp,material='Al')
    dP_damp = 1./(2.0*Qi*dQiInv_dNqp*Nqp/P)
    dP_dphs = 1./(-4.0*Qi/fr*dfr_dNqp*Nqp/P)

    return _calc_nep(psd,dP_damp,dP_dphs,tqp,tres)

def _calc_nep(psd, dP_damp, dP_dphs, tqp, tres):

    tres_term = np.sqrt(1 + (2*np.pi*psd.f*tres)**2)
    tqp_term  = np.sqrt(1 + (2*np.pi*psd.f*tqp)**2)

    tmpdict = {}
    tmpdict['frq']  = np.array(psd.f)
    tmpdict['amp']  = np.array(np.sqrt(psd.amplitude)*(dP_damp)*tres_term*tqp_term)
    tmpdict['phs']  = np.array(np.sqrt(psd.phase)*(dP_dphs)*tres_term*tqp_term)

    nepdata = PSDdata(tmpdict['frq'],tmpdict['amp'],tmpdict['phs'])

    return nepdata


