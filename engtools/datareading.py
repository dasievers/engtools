#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File Import Tools
"""

import pandas as pd
from io import StringIO
try:
    from timing import Timer
except ModuleNotFoundError:
    # Timer is helpful, but not necessary; create dummy class if not exist
    class Timer:
        def dummy(*args, **kwargs): pass
        def __getattr__(self, _): return self.dummy
import os


# =============================================================================
# Loading Data
# =============================================================================

def read_data(datapath, datetimecols=[0], tz=None, selectcols=None, dropcols=None,
                deprecated_names={},  overwrite_nulls=False, calcs=None):

    if '.xls' in datapath:
        importdf = pd.read_excel(datapath, parse_dates=datetimecols, index_col=0)
    else:
        with open(datapath, 'r') as f:
            s = f.read()
        if '\x00' in s:
            print('null characters detected in {}'.format(os.path.split(datapath)[1]))
            # delete "null" characters (from some Opto22 datalogs)
            s = s.replace('\x00', '')
            if overwrite_nulls:   #overwrite original file with nulls removed
                with open(datapath, 'w') as f:
                    f.write(s)
                print('original file replaced with clean version')
        if '\t' in s:
            delimiter = '\t'
        elif ',' in s:
            delimiter = ','
        else:
            delimiter = None  # error
        importdf = pd.read_csv(StringIO(s), delimiter=delimiter,
                               parse_dates=datetimecols, index_col=0)

    if tz is not None:
        # set the time zone
        importdf.index = importdf.index.tz_localize(tz)

    # rename and drop columns as specified
    importdf.rename(columns=deprecated_names, inplace=True)
    if selectcols is not None:
        dropcols = set(importdf.columns.values) - selectcols
    if dropcols is not None:
        importdf.drop(list(dropcols), axis=1, errors='ignore', inplace=True)


    # Additional Processing (if called for)
    if calcs is not None:
        importdf = calcs(importdf)

    return importdf


def batch_read(datapath_list, **kwargs):
    print('loading data files...')
    datalist = []
    timeit = Timer()
    for datapath in datapath_list:
        print(os.path.split(datapath)[1])
        importdf = read_data(datapath, **kwargs)
        datalist.append(importdf)

    masterdata = pd.concat(datalist, sort=True)
    masterdata.sort_index(inplace=True)


    if any(masterdata.index.duplicated()):  # need to collapse duplicate rows
        masterdata = masterdata.pivot_table(index=masterdata.index, dropna=False)  # this seems to work

    timeit.split('elapsed for import')
    return masterdata



