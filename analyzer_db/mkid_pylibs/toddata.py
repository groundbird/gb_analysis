################################################################
# TOD data
################################################################
from .basedata import BaseData
from .timearray import TimeArray
from .kiddata import KidData
from .psddata import PsdData
from . import misc

class TodData(BaseData, TimeArray, KidData):
    """
    Class for TOD data: TimeArray + KidData.
    time[sec] vs IQdata

    - **TodData(TimeArray,KidData,frq)**
    - **TodData(index,sps,kind,IQ-tuple,frq)**
    - **TodData(timestamp,kind,IQ-tuple,frq)**

    For TimeArray,
    :param sps: float of sampling rate [SPS(=Hz)]
    :param index: 1-D array of integer for timestamp
    :param timestamp: 1-D array of timestamps [sec]

    For KidData,
    :param IQ-tuple: array of measured data = (I,Q)
    :param kind: kind of data

    :param frq: float of signal frequency
    :param info: optional information
    """
    def __init__(self, *args, **kws):
        self._psd = None
        self._f = None

        if 'info' in kws:
            info = kws['info']
        else:
            info = dict()

        if not args:
            return

        self._f = args[-1]
        if len(args)==3:
            xdata = args[0] # TimeArray
            ydata = args[1] # KidData
        elif len(args)==4:
            xdata = TimeArray(args[0]) # timestamp
            ydata = KidData(args[1], args[2]) # kind + IQ
        elif len(args)==5:
            xdata = TimeArray(args[0],args[1]) # index + sps
            ydata = KidData(args[2], args[3]) # kind + IQ
        else:
            return

        super().__init__((xdata,ydata),info)

    @property
    def f(self):
        "Frequency [Hz]"
        return self._f
    @property
    def frequency(self):
        "Frequency [Hz]"
        return self._f
    @property
    def fGHz(self):
        "Frequency [GHz]"
        return self._f/1e9
    @property
    def fMHz(self):
        "Frequency [MHz]"
        return self._f/1e6

    @property
    def psd(self):
        "PSD"
        return self._psd

    def add_psddata(self, frq, amp, phs):
        """
        add PSD data
        :param frq: list of PSD frequency
        :param amp: list of PSD amplitude
        :param phs: list of PSD phase
        """
        self._psd = PsdData(frq,amp,phs)

    def f_Hz(self,unit="G"):
        "Frequency [Hz]"
        if unit[-2:].lower() == 'hz':
            unit = unit[:-2]
        return self._f/misc.convunit(unit)

    def add_offset(self, offset):
        """
        add offset to the frequency array
        :param offset: frequency offset [Hz]
        """
        if offset is None: return None
        self._f = self._f + misc.tofrq(offset)
        return self._f

