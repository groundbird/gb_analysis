import numpy as np

from . import misc

class PsdAnalyzer(object):
    def __init__(self, swp=None, cswp=None, tod=None, ctod=None, kidname='MKID'):
        '''
        Utilizer for KID SPD analysis
        :param Swpdata swp: Swpdata for on-resonance frq
        :param Swpdata cswp: Swpdata for off-resonance frq
        :param list[TodData] tod: 1D list of TodData for on-resonance frq
        :param list[TodData] ctod: 1D list of TodData for off-resonance frq
        :param str kidname: optional name for your use
        '''
        self._swpdata = swp
        self._cont_swpdata = cswp
        self._toddatalist = tod
        self._cont_toddatalist = ctod
        self._psd = None
        self._nep = None
        self._psdfitresult = None
        self._trgfitresult = None

    @property
    def swp(self):
        return self._swpdata

    @property
    def cswp(self):
        return self._cont_swpdata

    @property
    def tod(self):
        return self._toddatalist[0]

    @property
    def ctod(self):
        return self._cont_toddatalist[0]

    @property
    def todlist(self):
        return self._toddatalist

    @property
    def ctodlist(self):
        return self._cont_toddatalist

    @property
    def psd(self):
        return self._psd

    @property
    def psdfitresult(self):
        return self._psdfitresult

    @property
    def tqp(self):
        if self._psdfitresult is not None:
            return self.tqp_psd
        elif self._trigfitresult is not None:
            return self.tqp_trig
        else:
            return None

    @property
    def tqp_psd(self):
        if self._psdfitresult is None:
            return None
        return self._psdfitresult.params['tqp'].value

    @property
    def tqp_trig(self):
        if all(tod.fitresult is None for tod in self._toddatalist):
            return None
        ## This may have to be tuned.
        return np.array([tod.fitresult.params['tqp'].value for tod in self._todlist if tod.fitresult is not None]).mean()

    @property
    def nep(self):
        return self._nep

    def fitIQ(self, **kws):
        if self._swpdata is not None:
            self._swpdata.fitIQ(**kws)

        from .modIQ import rwIQ, modIQ
        from .kiddata import KidData
        if self._cont_toddatalist is not None:
            for tod,ctod in zip(self._toddatalist,self._cont_toddatalist):
                miq = modIQ(tod.iq,ctod.iq)
                tod.mdata = KidData('I-Q', (miq.real,miq.imag))
                if self._swpdata is not None:
                    riq = rwIQ(tod.f,tod.iq,self._swpdata.fitresult)
                    tod.rwdata = KidData('I-Q', (riq.real,riq.imag))
                    riq = rwIQ(ctod.f,ctod.iq,self._swpdata.fitresult)
                    ctod.rwdata = KidData('I-Q', (riq.real,riq.imag))
                    rmiq = rwIQ(tod.f,tod.mdata.iq,self._swpdata.fitresult)
                    tod.rwmdata = KidData('I-Q', (rmiq.real,rmiq.imag))
                pass
            pass
        else:
            for tod in self._toddatalist:
                if self._swpdata is not None:
                    riq = rwIQ(tod.f,tod.iq,self._swpdata.fitresult)
                    tod.rwdata = KidData('I-Q', (riq.real,riq.imag))
                pass
            pass

    ### Trigger ######

    def anatrig(self, ftype, inpath, frq, index=0, cont_index=None, kind=None, silent=True, keepall=False):
        '''
        Wrapper for analysis of events triggered by cosmic-ray muons.
        This provides tauqp via their signal-shape fittings.
        See fit_trigdata_readfromfile() in trig.py for option descriptions.
        '''
        from .trig import fit_trigdata_readfromfile
        self._trgfitresult = fit_trigdata_readfromfile(self._swpdata.fitresult, ftype, inpath, frq, index, cont_index, kind, silent, keepall)
        pass

    ### PSD ##########
    def calcpsd(self, kind=None, fold=None, dofit=True, silent_fit=True, is_norm=False, customfold=None,customlim=None,customind=None):
        '''
        Wrapper for PSD calculation.
        :param str kind: TOD data type to be calculated. One of ['rwmdata', 'rwdata', 'mdata', 'data'].
                         'rw': rewound one. 'm': modified by using off-resonance.
                         If not specified, the priority is 'rwmdata' --> 'rwdata' --> 'mdata' --> 'data'.
        :param fold: Folded-number for the PSD calculation. Lower number will provide more precise spectrum.
                     If fold == 'custom', the magic numbers for 3 tod data (rate == 1k, 100k, 1M SPS) will be used.
                     If fold is the list of int, different folded-number will be used for each data in self.tod.
        :param bool dofit: PSD roll-off fit will be performed, which provides tauqp.
        '''
        from .psd import calc_psd_fromdata,fit_psd
        from .psddata import PsdData
        rates = np.array([toddata.rate for toddata in self._toddatalist])
        f_lims = np.array([0.1] + (rates/2.).tolist())
        f_inds = np.arange(len(self._toddatalist))

        f_folds = fold
        if not hasattr(fold,'__iter__'):
            f_folds = [fold] * len(self._toddatalist)

        tmpdict={}
        tmpdict['frq'] = []
        tmpdict['amp'] = []
        tmpdict['phs'] = []

        if customfold is not None and customlim is not None and customind is not None:
            f_lims = customlim
            f_inds = customind
            f_folds = customfold
            for i,ind in enumerate(f_inds):
                f,a,p = calc_psd_fromdata(self._toddatalist[ind], kind=kind, fold=f_folds[i], is_norm=is_norm)
                sel = (f_lims[:-1][i] < f) & (f < f_lims[1:][i])
                #sel = misc.widener(sel)
                tmpdict['frq'] += f[sel].copy().tolist()
                tmpdict['amp'] += a[sel].copy().tolist()
                tmpdict['phs'] += p[sel].copy().tolist()
        else:
            for i,ind in enumerate(f_inds):
                f,a,p = calc_psd_fromdata(self._toddatalist[ind], kind=kind, fold=f_folds[i], is_norm=is_norm)
                self._toddatalist[ind].add_psddata(f,a,p)
                sel = (f_lims[:-1][i] < f) & (f < f_lims[1:][i])
                #sel = misc.widener(sel)
                tmpdict['frq'] += self._toddatalist[ind].psd.f[sel].copy().tolist()
                tmpdict['amp'] += self._toddatalist[ind].psd.amplitude[sel].copy().tolist()
                tmpdict['phs'] += self._toddatalist[ind].psd.phase[sel].copy().tolist()

        self._psd = PsdData(tmpdict['frq'],tmpdict['amp'],tmpdict['phs'])
        if dofit: self._psdfitresult = fit_psd(self._psd, swpfitresult=self._swpdata.fitresult, fullfit=True, silent=silent_fit)

    ### TempSweep ####
    def calcnep(self,tswp,V):
        fitres__dfr_dnqp = tswp.fit_linear('rolloff_Nqp','swp_fr',plotted=True,fitname='lin')[0]
        fitres__dQi_dnqp = tswp.fit_linear('rolloff_Nqp','swp_Qi',y_invert=True,plotted=True,fitname='lin')[0]
        dQiInv_dNqp = fitres__dQi_dnqp.params['a'].value
        dfr_dNqp = fitres__dfr_dnqp.params['a'].value

        Nqp = misc.convert_tqp_nqp(self.tqp)*V

        from .nep import calc_nep
        self._nep = calc_nep(self.psd, dfr_dNqp=dfr_dNqp, dQiInv_dNqp=dQiInv_dNqp, swpfitresult=self._swpdata.fitresult, tqp=self.tqp, Nqp=Nqp)
