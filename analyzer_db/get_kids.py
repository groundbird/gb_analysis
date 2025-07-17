#!/usr/bin/env python
'''Helper function to get KIDs data using DB information
'''
from pathlib import Path
import os
import sys
import django

from . import mkid_pylibs as klib
from .mkid_pylibs.readfile import readfile_swp, readfile_tod
from .kidslist import KidsList

MEASWEB_PATH = Path('/data/gb/db/meas_web')
sys.path.append(str(MEASWEB_PATH))

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "meas_web.settings")
django.setup()

from file_manager.models import Measurement

def get_kids(swp_mod, tod_mod, kidsfile):
    '''Obtain list of KIDs from Sweep object, Tod object and KidsList
    Parameters
    ----------
    swp_mod: file_manager.models.Sweep
        Sweep object
    tod_mod: file_manager.models.Tod
        Tod object
    kidsfile: kidslist.KidsList
        KidsList object

    Returns
    -------
    kids: list of klib.kidana.KidAnalyzer
        KIDs
    '''
    swps = readfile_swp('rhea',
                        filename=swp_mod.path,
                        index=-1, # all
                        lo=swp_mod.rheasweep.sg)

    tods = readfile_tod('rhea',
                        filename=tod_mod.path,
                        index=-1, # all
                        lo=tod_mod.rheatod.sg)

    kids_index = kidsfile.kids_index
    blinds_index = kidsfile.blinds_anaindex
    kids = []

    if len(kids_index) != len(blinds_index):
        raise Exception(f'Different length of indecise in `kids_index`({len(kids_index)}) and `blinds_index`({len(blinds_index)}). Currently one blind tone is supported for each KID.')

    for k_ind, b_ind in zip(kids_index, blinds_index):
        kid = klib.kidana.KidAnalyzer(swp=swps[k_ind],
                                      cswp=swps[b_ind],
                                      tod=tods[k_ind],
                                      ctod=tods[b_ind])
        kids.append(kid)

    return kids


def get_kids_id(meas_id):
    '''Obtain list of KIDs from measurement ID using the database
    Parameter
    ---------
    meas_id: int
        Measurement ID

    Returns
    -------
    kids: list of mkid_pylibs.kidana.KidAnalyzer
        KIDs
    '''
    meas_obj = Measurement.objects.filter(pk=meas_id).first()
    swp_objs = meas_obj.sweeps.all()
    tod_objs = meas_obj.tods.all()

    if swp_objs.count() != 1:
        print(f'# of sweeps: {swp_objs.count()}')
        print('Using the last one')

    swp_obj = swp_objs.last()

    if tod_objs.count() != 1:
        print(f'# of TODs: {tod_objs.count()}')
        print('Using the last one')

    tod_obj = tod_objs.last()

    if meas_obj.klpath != '':
        kidslist = KidsList(meas_obj.klpath)
    else:
        # Assume KidsList is at the same directory as tod
        tod_parent = Path(tod_obj.path).parent
        path_cand = tod_parent.glob('*.list')
        for _p in path_cand:
            _k = KidsList(_p)
            if _k.array_name == meas_obj.array:
                kidslist = _k
                break

    if kidslist is None:
        raise Exception('Kidslist not found')

    return get_kids(swp_obj, tod_obj, kidslist)
