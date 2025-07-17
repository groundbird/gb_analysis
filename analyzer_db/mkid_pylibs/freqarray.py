import copy
import numpy as np
from . import misc

################################################################
# freq data for indexing
# should have x (as property), __len__
################################################################

class FreqArray(object):
    """
    A class for representing array of frequency

    - **FreqArray(frqs)**:
    :param frqs: a 1-D array of frequency [Hz]
    """
    def __init__(self, frqs):
        self._data = np.asarray(frqs)
    def __len__(self):
        return len(self._data)
    def __getitem__(self, key):
        c = copy.copy(self)
        c._data = self._data[key]
        return c
    def fields(self):
        "data names inside object"
        return ('f',)
    def unpack(self):
        "return list of arrays in the same order with fields()"
        return [getattr(self, k) for k in self.fields()]

    @property
    def x(self):
        "Frequency [Hz]"
        return self._data
    @property
    def f(self):
        "Frequency [Hz]"
        return self._data
    @property
    def fGHz(self):
        "Frequency [GHz]"
        return self._data/1e9
    @property
    def fMHz(self):
        "Frequency [MHz]"
        return self._data/1e6
    @property
    def frequency(self):
        "Frequency [Hz]"
        return self._data

    def f_Hz(self,unit="G"):
        "Frequency [Hz]"
        if unit[-2:].lower() == 'hz':
            unit = unit[:-2]
        return self._data/misc.convunit(unit)

    def add_offset(self, offset):
        '''
        add offset to the frequency array
        :param offset: frequency offset [Hz]
        '''
        if offset is None: return None
        self._data = self._data + misc.tofrq(offset)
        return self._data

