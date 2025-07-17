################################################################
# Sweep data
################################################################
from .basedata import BaseData
from .freqarray import FreqArray
from .kiddata import KidData
from . import kidfit

class SwpData(BaseData, FreqArray, KidData):
    """
    Class for sweep data: FreqArray + KidData
    freq[Hz] vs IQdata

    - **SwpData(FreqArray,KidData)**:
    - **SwpData(frq,kind,IQ-tuple)**:
    - **SwpData(frq,offset,kind,IQ-tuple)**:
    :param frq: 1-D array of float for frequency [Hz]
    :param offset: float of offset frequency [Hz]
    :param IQ-tuple: array of measured data = (I,Q)
    """
    def __init__(self, *args, **kws):
        if 'info' in kws:
            info = kws['info']
        else:
            info = dict()
        if not args:
            return
        if len(args)==2:
            xdata = args[0] # FreqArray
            ydata = args[1] # KidData
        elif len(args)==3:
            xdata = FreqArray(args[0]) # frq
            ydata = KidData(args[1],args[2]) # kind + IQ
        elif len(args)==4:
            xdata = FreqArray(args[0]) # frq
            xdata.add_offset(args[1]) # offset
            ydata = KidData(args[2],args[3]) # kind + IQ
        else:
            return
        super().__init__((xdata, ydata), info)

    @property
    def fitdata(self):
        return self._fitdata
    @fitdata.setter
    def fitdata(self, val):
        self._fitdata = val

    def fitIQ(self, ciq=None, **kws):
        '''
        Fit registered swpdata

        :param ciq: off-resonance IQ
        :param kws: directly passed as `kidfit.fitIQ_onepeak(self, **kws)`
        '''
        self._fitresult = kidfit.fitIQ_onepeak(self, **kws)
        fiq = self._fitresult.fitted(self.f)
        self.fitdata = KidData('I-Q', (fiq.real,fiq.imag))
        from .modIQ import rwIQ, modIQ
        riq = rwIQ(self.f,self.iq,self._fitresult)
        self.rwdata = KidData('I-Q', (riq.real,riq.imag))
        rfiq = rwIQ(self.f,fiq,self._fitresult)
        self.fitdata.rwdata = KidData('I-Q', (rfiq.real,rfiq.imag))
        if ciq is not None:
            miq = modIQ(self.iq,ciq)
            self.mdata = KidData('I-Q', (miq.real,miq.imag))
            rmiq = rwIQ(self.f,self.mdata.iq,self._fitresult)
            self.rwmdata = KidData('I-Q', (rmiq.real,rmiq.imag))
            pass
        return self._fitresult

        #from .modIQ import add_moddata
        #add_moddata(self.f,self.fitdata,self.fitresult)
        #add_moddata(self.f,self.data,self.fitresult)

        # additional corrections by Kutsuma
        #self._fitresult.info['circle_fit'] = dict()
        #self._fitresult.info['circle_fit']['zcfit'], self._fitresult.info['circle_fit']['rfit'] = misc.circle_fit(self.rwdata.iq)
        #add_moddata(self.f,self.fitdata,self.fitresult)
        #add_moddata(self.f,self,self.fitresult)

