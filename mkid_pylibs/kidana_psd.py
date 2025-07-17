import numpy as np

from . import misc

class PSDAnalyzer(object):
    def __init__(self, swp=None, cswp=None, tod=None, ctod=None, kidname='MKID'):
        '''
        Utilizer for KID SPD analysis
        :param Swpdata swp: Swpdata for on-resonance frq
        :param Swpdata cswp: Swpdata for off-resonance frq
        :param list[TODdata] tod: 1D list of TODdata for on-resonance frq
        :param list[TODdata] ctod: 1D list of TODdata for off-resonance frq
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
        from .KidGlobalResponse import KidGlobalResponse
        if self._cont_toddatalist is not None:
            for tod,ctod in zip(self._toddatalist,self._cont_toddatalist):
                miq = modIQ(tod.iq,ctod.iq)
                tod.mdata = KidGlobalResponse('I-Q', (miq.real,miq.imag))
                if self._swpdata is not None:
                    riq = rwIQ(tod.f,tod.iq,self._swpdata.fitresult)
                    tod.rwdata = KidGlobalResponse('I-Q', (riq.real,riq.imag))
                    riq = rwIQ(ctod.f,ctod.iq,self._swpdata.fitresult)
                    ctod.rwdata = KidGlobalResponse('I-Q', (riq.real,riq.imag))
                    rmiq = rwIQ(tod.f,tod.mdata.iq,self._swpdata.fitresult)
                    tod.rwmdata = KidGlobalResponse('I-Q', (rmiq.real,rmiq.imag))
                pass
            pass
        else:
            for tod in self._toddatalist:
                if self._swpdata is not None:
                    riq = rwIQ(tod.f,tod.iq,self._swpdata.fitresult)
                    tod.rwdata = KidGlobalResponse('I-Q', (riq.real,riq.imag))
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
    def calcpsd(self, kind=None, fold=None, dofit=True, silent_fit=True):
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
        from .PSDdata import PSDdata
        rates = np.array([toddata.rate for toddata in self._toddatalist])

        if fold == "custom":
            ## magic numbers are available to smooth the PSD line! : assuming 3 psd data (1k, 100k, 1M SPS)
            f_folds = [32,64,32,64,32,64]
            f_lims = np.array([0.01, 2e2, 5e2, 2e4, 5e4, 2e5, 2e6])
        elif hasattr(fold,'__iter__'):
            ntod = len(self._toddatalist)
            if len(fold) >= ntod:
                f_folds = fold[:ntod]
            elif len(fold) < ntod:
                f_folds = list(fold) + [None]*(ntod-len(fold))
            f_lims = np.array([0.1] + (rates/2.).tolist())
        else:
            f_folds = [fold] * len(self._toddatalist)
            f_lims = np.array([0.1] + (rates/2.).tolist())

        f_ran = zip(f_lims[:-1], f_lims[1:])

        for tod,ifold in zip(self._toddatalist,f_folds):
            f,a,p = calc_psd_fromdata(tod, kind=kind, fold=ifold)
            tod.add_psddata(f,a,p)

        frqs  = [tod.psd.f for tod in self._toddatalist]
        amps  = [tod.psd.amplitude for tod in self._toddatalist]
        phss = [tod.psd.phase for tod in self._toddatalist]
        sel_inds = [(ran[0] < fs) & (fs < ran[1]) for fs, ran in zip(frqs, f_ran)]
        sel_inds = [misc.widener(sel_ind) for sel_ind in sel_inds]

        tmpdict={}
        tmpdict['frq'] = np.array(sum([frq[sel_ind].tolist() for frq, sel_ind in zip(frqs, sel_inds)], []))[1:]
        tmpdict['amp'] = np.array(sum([amp[sel_ind].tolist() for amp, sel_ind in zip(amps, sel_inds)], []))[1:]
        tmpdict['phs'] = np.array(sum([phs[sel_ind].tolist() for phs, sel_ind in zip(phss, sel_inds)], []))[1:]

        self._psd = PSDdata(tmpdict['frq'],tmpdict['amp'],tmpdict['phs'])

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

