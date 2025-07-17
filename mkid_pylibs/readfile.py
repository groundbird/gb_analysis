import warnings
from .TODdata import TODdata
from .Swpdata import Swpdata
from . import misc
from .rhea_comm.lib_read_rhea import read_rhea_swp, read_rhea_tod, read_rhea_tod_sync
from .parsefile import parse_vna_csv, parse_vna_citi,parse_riken_csv

def _readfile_swp_rhea(filename,index, IQkind = 'I-Q', **kws):
    tmpdata = read_rhea_swp(filename,ismult=True,miniret=True, **kws)
    if index == 'all': index = range(len(tmpdata))
    if IQkind == 'I-Q':
        print('IQkind = I-Q')
        return [Swpdata(tmpdata[ii]['freq'],'I-Q',(tmpdata[ii]['I'], tmpdata[ii]['Q'])) for ii in index]
    elif IQkind == 'I-Qerr':
        print('IQkind = I-Qerr')
        return [Swpdata(tmpdata[ii]['freq'],'I-Qerr',(tmpdata[ii]['I'], tmpdata[ii]['Q'], tmpdata[ii]['Ierr'], tmpdata[ii]['Qerr'])) for ii in index]
    else:
        print('Kind error, you should set kind is `I-Q` or `I-Qerr`')

def _readfile_tod_rhea(filename,index,sync,**kws):
    if sync:
        tmpdata = read_rhea_tod_sync(filename,miniret=True,**kws)
    else:
        tmpdata = read_rhea_tod(filename,miniret=True,**kws)
    if index == 'all': index = range(len(tmpdata))
    if sync:
        return [TODdata(tmpdata[ii]['time'], 'I-Q', (tmpdata[ii]['I'], tmpdata[ii]['Q']), tmpdata[ii]['freq'], info={'n_rot':tmpdata[ii]['n_rot'],'sync_off':tmpdata[ii]['sync_off']}) for ii in index]
    else:
        return [TODdata(tmpdata[ii]['time'], 'I-Q', (tmpdata[ii]['I'], tmpdata[ii]['Q']), tmpdata[ii]['freq']) for ii in index]

def _readfile_vna(filename,**kws):
    warnings.warn(f'This function is now outdated. Contact to shonda.')
    if filename[-4:] == '.csv':
        tmpinfo, tmpdata = parse_vna_csv(filename)
    elif filename[-4:] == '*.cti':
        tmpinfo, tmpdata = parse_vna_citi(filename)
    else:
        raise misc.MKIDDataException(f'Invalid file format: {filename[-4:]}')
    return tmpinfo,tmpdata

def _readfile_tod_vna(filename,index,**kws):
    tmpinfo, tmpdata = _readfile_vna(filename)
    if index == 'all': index = range((len(tmpdata)-1)//2)
    return [TODdata(tmpdata[0], 'dB-DEG', (tmpdata[ii*2+1], tmpdata[ii*2+2]), tmpinfo['CW Freq']) for ii in index]

def _readfile_swp_vna(filename,index,**kws):
    tmpinfo, tmpdata = _readfile_vna(filename)
    if index == 'all': index = range((len(tmpdata)-1)//2)
    return [Swpdata(tmpdata[0], 'dB-DEG', (tmpdata[ii*2+1], tmpdata[ii*2+2])) for ii in index]

def _readfile_tod_riken(filename,index,**kws):
    warnings.warn(f'This function is now outdated. Contact to shonda.')
    tmpinfo,tmpdata = parse_riken_csv(filename)
    if index == 'all': index = range((len(tmpdata)-1)//2)
    return [TODdata(tmpdata[0], 'dB-DEG', (tmpdata[ii*2+1], tmpdata[ii*2+2]), 0) for ii in index]

def _readfile_swp_riken(filename,index,**kws):
    warnings.warn(f'This function is now outdated. Contact to shonda.')
    tmpinfo,tmpdata = parse_riken_csv(filename)
    if index == 'all': index = range((len(tmpdata)-1)//2)
    return [Swpdata(tmpdata[0], 'dB-DEG', (tmpdata[ii*2+1], tmpdata[ii*2+2])) for ii in index]

def readfile_tod(kind, filename, index=0, lo=0, sync=True, **kws):
    """
    read TODdata from given file

    :param kind: one of ['rhea', 'riken', 'vna']
    :param filename: name of input file
    :param index: KID index number to be assigned in the input file
    """
    _todkinds = dict()
    _todkinds['riken'] = _readfile_tod_riken
    _todkinds['vna']   = _readfile_tod_vna
    _todkinds['rhea']  = _readfile_tod_rhea

    return _readfile(_todkinds,kind,filename,index,lo,sync=sync,**kws)

def readfile_swp(kind, filename, index=0, lo=0, **kws):
    """
    read Swpdata from given file

    :param kind: one of ['rhea', 'riken', 'vna']
    :param filename: name of input file
    :param index: KID index number to be assigned in the input file
    """
    _swpkinds = dict()
    _swpkinds['riken'] = _readfile_swp_riken
    _swpkinds['vna']   = _readfile_swp_vna
    _swpkinds['rhea']  = _readfile_swp_rhea

    return _readfile(_swpkinds,kind,filename,index,lo,**kws)

def _readfile(dictkinds, kind, filename, index, lo, **kws):
    is_iterative = True
    if not hasattr(index, '__iter__'):
        if index < 0:
            index = 'all'
        else:
            is_iterative = False
            index = [index]
            pass
        pass

    if not kind.lower() in dictkinds:
        raise misc.MKIDDataException(f'Invalid option `{kind}`. Choose from {dictkinds.keys()}.')

    ret = dictkinds[kind.lower()](filename,index, **kws)
    for rr in ret:
        rr.add_offset(lo)
        rr.info['filename'] = filename
        pass
    if not is_iterative: ret = ret[0]

    return ret

