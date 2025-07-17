import numpy as np
import pandas as pd
import datetime


def read_vna(path):
    data = pd.read_table(path, skiprows=1,sep=' ', header=None,names=['freq', 'I', 'Q', 'logmag', 'phase'], skipfooter=1, engine='python')
    date = datetime.datetime.strptime(path.split('.')[0].split('/')[-1], '%Y-%m%d-%H%M%S')
    return data, date


def read_mulvna(paths):
    datas = []
    dates = []
    for ipath in paths:
        data, date = read_vna(ipath)
        datas.append(data)
        dates.append(date)
    return datas, dates