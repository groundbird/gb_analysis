from functools import reduce
import numpy as np
import copy
import operator
from .Drawer import Drawer

###########################################################################
# class for combined IQ(amp/phs) data and indexing array (freq/time)
###########################################################################
class BaseMultiData(Drawer):
    """
    base class for combined data array

    datalist is a list of object with
    - __len__
    - __getitem__
    """
    def __init__(self, datalist, info=None, *args, **kws):
        if len(set(map(len,datalist))) != 1:
            print(datalist)
            raise RuntimeError('ERROR:: Given arrays are not of equal length')
        self._datalist = datalist
        if info is None:
            self.info = dict()
        else:
            self.info = info.copy()
        self._debug = False
        if 'debug' in kws:
            self._debug = True
    def __len__(self):
        return len(self._datalist)
    def __getitem__(self, key):
        c = copy.copy(self)
        c._datalist = [d[key] for d in self._datalist]
        return c
    def __getattr__(self, attr):
        ary_attr = attr.split('.')
        if ary_attr[0]  == '_datalist':
            raise AttributeError()
        tmpary = [getattr(d,attr) for d in self._datalist if hasattr(d,attr)]
        if len(tmpary) == 0:
            raise AttributeError(f'No candidates of {attr} in {self.__class__.__name__} ({self._datalist})')
#        elif len(tmpary) > 1 and self._dedug:
#            import warnings
#            warnings.warn(f'Multiple arrays for {attr} in {self.__class__.__name__} ({self._datalist}).')
        return tmpary[0]
    def fields(self):
        "data names inside object"
        return (d.fields() for d in self._datalist)
    def unpack(self):
        "return list of arrays in the same order with fields()"
        return [getattr(self, k) for k in self.fields()]

    def add_info(self, keyword, valueword):
        self.info[keyword] = valueword

    def print_info(self):
        for i, par in enumerate(self.info):
            print(f'{i:02d} : {par:15s} -> {self.info[par]}')

