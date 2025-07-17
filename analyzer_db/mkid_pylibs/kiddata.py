import numpy as np
from . import misc

## Data Definitions for Different Formats
class _KidGlobalResponse_IQData:
    """
    given normal IQ data
    :param data: array(0:I, 1:Q)
    """
    def __init__(self, data):
        self._data = dict()
        self._data['I'] = np.asarray(data[0])
        self._data['Q'] = np.asarray(data[1])
    def __len__(self):
        return len(self._data['I'])
    def __getitem__(self, key):
        return KidGlobalResponse('I-Q',(np.asarray(self._data['I'][key]),np.asarray(self._data['Q'][key])))
    def i(self):
        return self._data['I']
    def q(self):
        return self._data['Q']
    def iq(self):
        return self.i() + 1j*self.q()
    def deg(self):
        return np.arctan2(self.q(), self.i()) * 180/np.pi
    def amplitude(self):
        return np.sqrt(self.i()**2 + self.q()**2)
    def db(self):
        return misc.amplitude_to_dB(self.amplitude())
    def fields(self):
        return ('i', 'q')

class _KidGlobalResponse_IQGainData:
    """
    given IQdata with gain
    :param data: array(0:I, 1:Q, 2:gain)
    """
    def __init__(self, data):
        self._data = dict()
        self._data['I'] = np.asarray(data[0])
        self._data['Q'] = np.asarray(data[1])
        self._data['gain'] = data[2]
    def __len__(self):
        return len(self._data['I'])
    def __getitem__(self, key):
        return KidGlobalResponse('I-Q-Gain',(self._data['I'][key],self._data['Q'][key],self._data['gain'][key]))
    def i(self):
        return self._data['I'] / self._data['gain']
    def q(self):
        return self._data['Q'] / self._data['gain']
    def iq(self):
        return self.i() + 1j*self.q()
    def deg(self):
        return np.arctan2(self.q(), self.i()) * 180/np.pi
    def amplitude(self):
        return np.sqrt(self.i()**2 + self.q()**2)
    def db(self):
        return misc.amplitude_to_dB(self.amplitude())
    def fields(self):
        return ('i', 'q')

class _KidGlobalResponse_dBDEGData:
    """
    given dB and deg data
    :param data: array(0:dB, 1:deg)
    """
    def __init__(self, data):
        self._data = dict()
        self._data['dB']  = np.asarray(data[0])
        self._data['DEG'] = np.asarray(data[1])
    def __len__(self):
        return len(self._data['dB'])
    def __getitem__(self, key):
        return KidGlobalResponse('dB-DEG',(self._data['dB'][key],self._data['DEG'][key]))
    def i(self):
        if not 'I' in self._data:
            self._data['I'] = misc.dB_to_amplitude(self.db())*np.cos(self.deg()/180.0*np.pi)
        return self._data['I']
    def q(self):
        if not 'Q' in self._data:
            self._data['Q'] = misc.dB_to_amplitude(self.db())*np.sin(self.deg()/180.0*np.pi)
        return self._data['Q']
    def deg(self):
        return self._data['DEG']
    def iq(self):
        return self.i() + 1j*self.q()
    def amplitude(self):
        return misc.dB_to_amplitude(self.db())
    def db(self):
        return self._data['dB']
    def fields(self):
        return ('db', 'deg')

class _KidGlobalResponse_dBData():
    """
    given only dB data (unavailable phase)
    :param data: array(0:dB)
    """
    def __init__(self, data):
        self._data = dict()
        self._data['dB']  = np.asarray(data[0])
    def __len__(self):
        return len(self._data['dB'])
    def __getitem__(self, key):
        return KidGlobalResponse('dB',(self._data['dB'][key]))
    def i(self):
        raise misc.MKIDDataException('no phase data!')
    def q(self):
        raise misc.MKIDDataException('no phase data!')
    def iq(self):
        raise misc.MKIDDataException('no phase data!')
    def deg(self):
        raise misc.MKIDDataException('no phase data!')
    def amplitude(self):
        return misc.dB_to_amplitude(self.db())
    def db(self):
        return self._data['dB']
    def fields(self):
        return ('db',)

_kiddatakinds = dict()
_kiddatakinds['I-Q']       = _KidGlobalResponse_IQData
_kiddatakinds['dB-DEG']    = _KidGlobalResponse_dBDEGData
_kiddatakinds['I-Q-Gain']  = _KidGlobalResponse_IQGainData
_kiddatakinds['dB']        = _KidGlobalResponse_dBData

class KidData(object):
    """
    A class to store IQdata providing vairables:
    `i`, `q`, `iq`, `amplitude`, `db`, `deg`, `phase`

    - **KidGlobalResponse(kind, filename)**:
    :param kind: one of ['I-Q', 'dB-DEG', 'I-Q-Gain', 'dB']
    :param data: tuple of data, corresponding to kind
    """
    def __init__(self, kind, data):
        self._kind      = kind
        self._rwdata    = None
        self._mdata     = None
        self._rwmdata   = None
        self._fitresult = None
        #if len(set(map(len,data))) != 1:
        #   raise RuntimeError('ERROR:: Given arrays are not of equal length')
        if not kind in _kiddatakinds:
            raise misc.MKIDDataException('data format not implemented')
        self._data = _kiddatakinds[kind](data)
    def __len__(self):
        return len(self._data)
    def __getitem__(self, key):
        return self._data[key]
    def fields(self):
        "data names inside object"
        return self._data.fields()
    def unpack(self):
        "return list of arrays in the same order with fields()"
        return [getattr(self, k) for k in self.fields()]
    def xrange(self, beg, end):
        return np.where((beg <= self.x) and (self.x <= end))

    @property
    def i(self):
        "I"
        return self._data.i()
    @property
    def q(self):
        "Q"
        return self._data.q()
    @property
    def iq(self):
        "IQ"
        return self._data.iq()
    @property
    def amplitude(self):
        "S21 in amplitude"
        return self._data.amplitude()
    @property
    def coramplitude(self):
        "normalized S21 in amplitude"
        return self._data.amplitude()*2
    @property
    def db(self):
        "S21 [dB]"
        return self._data.db()
    @property
    def deg(self):
        "S21 in phase [deg]"
        return self._data.deg()
    @property
    def phase(self):
        "S21 in phase [rad:(-pi,pi)]"
        return np.angle(self._data.iq())
    @property
    def corphase(self):
        "S21 in phase [rad]"
        angletmp = -1*np.angle(self._data.iq())
        #if np.any(np.abs(angletmp[:-1]-angletmp[1:])>2.8):
        angletmp = angletmp - np.sign(angletmp)*np.pi
        return angletmp
    @property
    def lincorphase(self):
        "linearlized phase"
        return misc.linearized_phase(self.corphase)

    @property
    def data(self):
        "IQdata"
        return self
    @property
    def rwdata(self):
        "Rewound IQdata"
        return self._rwdata
    @rwdata.setter
    def rwdata(self,val):
        self._rwdata = val
    @property
    def mdata(self):
        "Modified IQdata subtracted by blind tone"
        return self._mdata
    @mdata.setter
    def mdata(self,val):
        self._mdata = val
    @property
    def rwmdata(self):
        "Rewound+modified IQdata"
        return self._rwmdata
    @rwmdata.setter
    def rwmdata(self,val):
        self._rwmdata = val

    @property
    def fitresult(self):
        return self._fitresult

    def add_IQdata(self,iq,kind):
        """
        Add IdQdata
        :param iq: IQdata list
        :param kind: the member name to be defined (`data`, `rwdata`, `mdata`, `rwmdata`)
        """
        setattr(self,'_'+kind,KidGlobalResponse('I-Q', (iq.real,iq.imag)))
        pass

