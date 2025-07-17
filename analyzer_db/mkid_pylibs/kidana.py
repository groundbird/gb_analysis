
class KidAnalyzer(object):
    def __init__(self, swp=None, cswp=None, tod=None, ctod=None, kidname='MKID'):
        '''
        Utilizer for KID analysis
        :param Swpdata swp: Swpdata for on-resonance frq
        :param Swpdata cswp: Swpdata for off-resonance frq
        :param list[TodData] tod: 1D list of TodData for on-resonance frq
        :param list[TodData] tod: 1D list of TodData for off-resonance frq
        :param str kidname: optional name for your use
        '''
        self._swpdata = swp
        self._cont_swpdata = cswp
        self._toddata = tod
        self._cont_toddata= ctod
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
        return self._toddata

    @property
    def ctod(self):
        return self._cont_toddata

    @property
    def psd(self):
        return self._toddata.psd

    @property
    def psdfitresult(self):
        return self._psdfitresult

    @property
    def trgfitresult(self):
        return self._trgfitresult

    def fitIQ(self, **kws):
        from .modIQ import rwIQ, modIQ
        from .kiddata import KidData

        # on/off subtruction
        if self._cont_swpdata is not None:
            miq = modIQ(self._swpdata.iq,self._cont_swpdata.iq)
            self._swpdata.mdata = KidData('I-Q',(miq.real,miq.imag))

        self._swpdata.fitIQ(**kws) # fit with rawdata

        if self._swpdata.mdata is not None:
            rmiq = rwIQ(self._swpdata.f,self._swpdata.mdata.iq,self._swpdata.fitresult)
            self._swpdata.rwmdata = KidData('I-Q', (rmiq.real,rmiq.imag))

        if self._toddata is not None:
            riq = rwIQ(self._toddata.f,self._toddata.iq,self._swpdata.fitresult)
            self._toddata.rwdata = KidData('I-Q', (riq.real,riq.imag))

            if self._cont_toddata is not None:
                criq = rwIQ(self._cont_toddata.f,self._cont_toddata.iq,self._swpdata.fitresult)
                self._cont_toddata.rwdata = KidData('I-Q', (criq.real,criq.imag))

                miq = modIQ(self._toddata.iq,self._cont_toddata.iq)
                self._toddata.mdata = KidData('I-Q', (miq.real,miq.imag))
                rmiq = rwIQ(self._toddata.f,self._toddata.mdata.iq,self._swpdata.fitresult)
                self._toddata.rwmdata = KidData('I-Q', (rmiq.real,rmiq.imag))

    ### PSD ##########

    def calcpsd(self, kind=None, fold=None, dofit=True,varname='phase',fitsilent=False):
        '''
        Wrapper for PSD calculation.
        :param str kind: TOD data type to be calculated. One of ['rwmdata', 'rwdata', 'mdata', 'data'].
                         'rw': rewound one. 'm': modified by using off-resonance.
                         If not specified, the priority is 'rwmdata' --> 'rwdata' --> 'mdata' --> 'data'.
        :param fold: Folded-number for the PSD calculation. Lower number will provide more precise spectrum.
                     If fold == 'custom', the magic numbers for 3 tod data (rate == 1k, 100k, 1M SPS) will be used.
                     If fold is the list of int, different folded-number will be used for each data in self.tod.
        :param bool dofit: PSD roll-off fit will be performed, which provides tauqp.
        :apram str varname: 'phase' or 'amplitude' for PSD fitting
        '''
        from .psd import calc_psd_fromdata,fit_psd

        f,a,p = calc_psd_fromdata(self._toddata, kind=kind, fold=fold)
        self._toddata.add_psddata(f,a,p)
        if dofit: self._psdfitresult = fit_psd(self._toddata.psd, swpfitresult=self._swpdata.fitresult, fullfit=True, varname=varname,silent=fitsilent)

    ### Trigger ######

    def fittrig(self, kind=None, silent=True):
        '''
        Wrapper for analysis of events triggered by cosmic-ray muons.
        This provides tauqp via their signal-shape fittings.
        See fit_trigdata_readfromfile() in trig.py for option descriptions.
        '''
        from .trig import fit_trig
        self._trgfitresult = fit_trig(self._toddata, kind=kind, silent=silent)
        pass

    ### TempSweep ####

    # def anasweep(self, ftype, inpath_swp, frq, index=0, cont_index=None, inpath_tod=None, temps=None, V=None):
    #     '''
    #     Wrapper for temperature sweep measurements.
    #     :param str ftype: one of ['rhea', 'riken', 'vna']
    #     :param list[str] inpath_swp: filename for Swpdata (one filename for each temperature)
    #     :param list[str] inpath_tod: filename for TodData (one filename for each temperature)
    #     :param int index: index number for on-resonance
    #     :param int cont_index: index number for off-resonance
    #     :param list[float] temps: 1D list of temperature
    #     :param float V: KID volume [um3]
    #     '''
    #     from .sweep_meas import SweepMeas
    #     tswp = SweepMeas()
    #     if temps is not None:
    #         tswp.sweepVar(temps, 'temp')
    #         tswp.sweepVar(misc.convert_Kelvin_nqp(np.array(temps)), 'nqp')
    #         if V is not None:
    #             tswp.sweepVar(misc.convert_Kelvin_nqp(np.array(temps)) * V, 'Nqp')
    #     tswp.sweepVar_readfromfile('rhea',inpath_swp,4.99e9,0,1,inpath_tod,'psd_phs',label='temp[K]',nmax=1e5, plotted=True)
    #     tswp.sweepVar(misc.convert_tqp_nqp(tswp.get_vals('psd_phs_tqp')), 'rolloff_nqp')
    #     tswp.sweepVar(misc.convert_tqp_nqp(tswp.get_vals('psd_phs_tqp'))*V, 'rolloff_Nqp')
    #     return tswp

    # def calcnep(self,tswp,V):
    #     fitres__dfr_dnqp = tswp.fit_linear('rolloff_Nqp','swp_fr',plotted=True,fitname='lin')[0]
    #     fitres__dQi_dnqp = tswp.fit_linear('rolloff_Nqp','swp_Qi',y_invert=True,plotted=True,fitname='lin')[0]
    #     dQiInv_dNqp = fitres__dQi_dnqp.params['a'].value
    #     dfr_dNqp = fitres__dfr_dnqp.params['a'].value

    #     Nqp = misc.convert_tqp_nqp(self.tqp)*V

    #     from .nep import calc_nep
    #     self._nep = calc_nep(self.psd, dfr_dNqp=dfr_dNqp, dQiInv_dNqp=dQiInv_dNqp, swpfitresult=self._swpdata.fitresult, tqp=self.tqp, Nqp=Nqp)

