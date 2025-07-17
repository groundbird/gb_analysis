#!/usr/bin/env python3
'''Utility for analyzer_db
'''
from datetime import datetime, timezone, timedelta
from pathlib import Path
import sys
import os
import django

from .rhea_comm.lib_read_rhea import TodHeader

# Django preparation
MEASWEB_PATH = Path('/data/gb/db/meas_web')
sys.path.append(str(MEASWEB_PATH))
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "meas_web.settings")
django.setup()

from file_manager.models import Tod

def todpath2span(path):
    '''Get infromation from DB or TodHeader to get start and end datetime of the measurement
    Parameters
    ----------
    path: pathlib.Path
        TOD file path

    Returns
    -------
    dt_st: datetime.datetime
    '''
    dt_en = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc)

    tod_objs = Tod.objects.filter(path=f'{path}')
    if len(tod_objs) == 1:
        return tod_objs[0].mtime, dt_en

    # Reaches here if no info in DB
    dur = TodHeader(path).duration
    dt_st = dt_en - timedelta(seconds=dur)
    return dt_st, dt_en
