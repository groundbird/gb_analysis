from .Drawer import Drawer

################################################################
# general array class with string format
################################################################

_general_array_template = \
"""
import numpy as np
import copy

class {typename}(Drawer):
    '''
    array class of {fields}.
    '''
    def __init__(self, arr, info=None):
        '''
        :param arr: a n-D array
        :param info: (optional) a dict for attaching information
        '''
        if np.ndim(arr)==1:
            arr = [arr]
        if len(set(map(len,arr))) != 1:
            raise RuntimeError('ERROR:: Given arrays are not of equal length')
        self._datagen = np.asarray(arr)
        if info is None:
            self.info = dict()
        else:
            self.info = info.copy()
    def __len__(self):
        return len(self._datagen[0])
    def __getitem__(self, key):
        c = copy.copy(self)
        c._datagen = np.asarray(self._datagen[:, key])
        return c
    def fields(self):
        return ({field_strs})
    def unpack(self):
        return [getattr(self, k) for k in self.fields()]
{field_defs}
"""
_field_template = \
"""
    @property
    def {name}(self):
        "{doc}"
        return self._datagen[{index}]
"""

def named_array(typename, field_names, docs=None, verbose=False):
    """
    general n-D array class

    :param typename: name for generated class
    :param field_name: a list of string. name for fields
    :param docs: a list of string. brief descriptions for fields. If not given, field_name is used.
    :param verbose: print class definition string just before it is actually defined.
    """
    if docs is None:
        field_defs = [_field_template.format(index=i,name=n, doc=n)
                      for i, n in enumerate(field_names)]
    else:
        field_defs = [_field_template.format(index=i,name=n, doc=d)
                      for i, (n, d) in enumerate(zip(field_names, docs))]

    field_strs = '("' + ('", "'.join(field_names)) + '",)'

    s = _general_array_template.format(typename=typename,
                                       fields=', '.join(field_names),
                                       field_defs='\n'.join(field_defs),
                                       field_strs=(field_strs))

    namespace = dict(__name__=f'namedarray_{typename}')
    namespace['Drawer'] = Drawer
    exec(s, namespace)
    result = namespace[typename]
    result._source = s
    if verbose:
        print((result._source))

    try: # for pickling to work: copied from namedtuple definition
        import sys as _sys
        result.__module__ = _sys._getframe(1).f_globals.get('__name__', '__main__')
    except (AttributeError, ValueError):
        pass

    return result

