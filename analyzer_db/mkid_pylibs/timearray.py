import datetime
import numpy as np
from . import misc

################################################################
# time data for indexing
# should have x (as property), __len__
################################################################

class TimeArray(object):
    """
    A class for representing array of time

    - **TimeArray(timestamp)**: time is given by timestamp itself.
    :param timestamp: a 1-D array of time in second.

    - **TimeArray(time, sps)**: time is given by index * samplerate.
    :param index: a 1-D array of timestamp index.
    :param sps: a 1-D array of sampling rate [SPS(=Hz)].
    """
    def __init__(self, *args):
        """
        Initialize either with
        - 1 arg: timestamp
        - 2 args: index and sps
        """
        if len(args) == 1:
            self._data = _TimeArray_Timestamp(*args)
        elif len(args) == 2:
            self._data = _TimeArray_Index(*args)
        else:
            raise misc.MKIDDataException('invalid initialization of TimeArray')
        pass
    def __len__(self):
        return len(self._data)
    def __getitem__(self, key):
        return self._data[key]
    def fields(self):
        "data names inside object"
        return ('t',)
    def unpack(self):
        "return list of arrays in the same order with fields()"
        return [getattr(self, k) for k in self.fields()]
    def xrange(self, beg, end):
        return np.where((beg <= self.x) and (self.x <= end))

    @property
    def x(self):
        "Time [s]"
        return self._data.value()
    @property
    def t(self):
        "Time [s]"
        return self._data.value()
    @property
    def time(self):
        "Time [s]"
        return self._data.value()
    @property
    def index(self):
        "Index Number"
        return self._data.index()
    @property
    def rate(self):
        "Rate [SPS]"
        return self._data.rate()

class _TimeArray_Index(object):
    def __init__(self, index, sps):
        self._index = np.asarray(index) #[non-units]
        self._rate = sps #[Hz]
    def __len__(self):
        return len(self._index)
    def __getitem__(self, key):
        return TimeArray(self._index[key], self._rate)
    def value(self):
        return self._index / self._rate
    def index(self):
        return self._index
    def rate(self):
        return self._rate

class _TimeArray_Timestamp(object):
    def __init__(self, timestamp):
        if type(timestamp[0]) is datetime.datetime:
            self._timestamp = np.asarray([(tt - timestamp[0]).total_seconds() for tt in timestamp])
        else:
            self._timestamp = np.asarray(timestamp) #[sec]
    def __len__(self):
        return len(self._timestamp)
    def __getitem__(self, key):
        return TimeArray(self._timestamp[key])
    def value(self):
        return self._timestamp
    def index(self):
        return np.arange(len(self._timestamp))
    def rate(self):
        if len(self._timestamp)>1:
            return 1./np.average(self._timestamp[1:] - self._timestamp[:-1]) #[Hz]

