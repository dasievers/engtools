#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Database management tools
"""

import pandas as pd
import os
import glob
import numpy as np
import sqlite3
from sqlalchemy import create_engine
from .datareading import batch_read

try:
    from .timing import Timer
except ModuleNotFoundError:
    # Timer is helpful, but not necessary; create dummy class if not exist
    class Timer:
        def dummy(*args, **kwargs): pass
        def __getattr__(self, _): return self.dummy
from warnings import warn


# =============================================================================
# Local Databases (sqlite, mysql, etc)
# =============================================================================

def dbupdate_sqlite(dbpath, datadir, table='values', **kwargs):
    """
    Update the sqlite database file with new data or create new if none exists.
        
    Parameters
    ----------
    dbpath : str
        Database file path. If none exists yet, define new path/name.
            
    datadir : str or list
        str: Path of the directory housing the files to be imported. Caution, all 
        files will be scanned so only put compatible files in this directory.
        list: Discrete list of complete path strings.
   
    Optional Parameters
    -------------------
    table : str
        Name of table in database to insert data. Defaults to `values`.
        
    **kwargs
        These kwargs are passed to read_data(). See more documentation there.
    """
#    presets = _presets(db, datadir)
    if type(datadir) == str:
        datapath_list = glob.glob(os.path.join(datadir, '*'))
    elif type(datadir) == list:
        datapath_list = datadir
    else:
        raise Exception('invalid datadir dtype')
    
    # check for new data
    # data will only be imported if not in the database already
    # 'import_history' table contains imported files and their date/time ranges
    
    # import previous file import history
    if os.path.isfile(dbpath):
        with sqlite3.connect(dbpath) as conn:
            c = conn.cursor()
            c.execute('SELECT "import_history"."file" FROM "import_history"')
            history = set(np.asarray(c.fetchall()).flatten())
        # delete entries from datapath_list (globbed) if also exist in history
        datapath_list = np.asarray(datapath_list)
        # print(datapath_list)
        newmask = np.ones(datapath_list.shape, dtype=bool)
        for i, f in enumerate(datapath_list):
            if os.path.split(f)[1] in history:
                newmask[i] = False
        datapath_list_new = datapath_list[newmask]
        # print(datapath_list)
    else:
        datapath_list_new = datapath_list
        
    # determine what columns already exist
    if os.path.isfile(dbpath):
        with sqlite3.connect(dbpath) as conn:
            c = conn.cursor()
            c.execute('SELECT name FROM sqlite_master WHERE type="table"')
            tables = set([t[0] for t in c.fetchall()])
            # can only check for colums if table exists
            if table in tables:
                c = conn.cursor()            
                c.execute('PRAGMA table_info("{}")'.format(table))
                exist_cols = {i[1] for i in c.fetchall()}
            else:
                exist_cols = set()
    else:
        exist_cols = set()

    # Import New Data
    if len(datapath_list_new) == 0:
        print('*** no new files were found to import ***')
        return
    
    timeit = Timer()    
    masterdata = batch_read(datapath_list_new, **kwargs)
    # change index to unix time integers (s) (for storage and query speed)
    masterdata.index = (masterdata.index.values.astype(float)/1e9).astype(int)
    
    # make a list of imported files to add to the db history table
    filelist = []
    for file in datapath_list_new:
        filelist.append(os.path.split(file)[1])
    filedata = pd.DataFrame(filelist, columns=['file'])
        
    timeit.split(report=False)        
    
    # Export to SQL
    # check for new colums and add them if needed (if table already exists)
    new_cols = list(set(masterdata.columns.values) - exist_cols)
    
    if os.path.isfile(dbpath):
        with sqlite3.connect(dbpath) as conn:
            c = conn.cursor()            
            c.execute('SELECT name FROM sqlite_master WHERE type="table"')
            tables = set([t[0] for t in c.fetchall()])
            # only modify if table exists
            if table in tables:
                for n in new_cols:
                    if 'int' in str(masterdata[n].dtype):
                        dtp = 'INT'
                    elif 'float' in str(masterdata[n].dtype):
                        dtp = 'FLOAT'
                    elif 'str' in str(masterdata[n].dtype):
                        dtp = 'VARCHAR(255)'
                    else:
                        dtp = ''
                    c.execute('ALTER TABLE "{0}" ADD "{1}" {2}'.format(
                            table, n, dtp))
           
    #TODO add functionality for using MySQL
    #engine = create_engine('mysql+pymysql://root:root@localhost:8889/'+switch)
    engine = create_engine('sqlite:///'+dbpath)
    masterdata.to_sql(table, engine, if_exists='append', index=True, index_label='date_time')
    filedata.to_sql("import_history", engine, if_exists='append', index=False)
    
    timeit.split('elapsed for SQL insert')



def query_sqlite(dbpath, table='values', start=None, end=None, tags=None, condition=None, addquery=None, 
                 tz='US/Mountain', verbose=False):
    """
    Query data from a SQLite database.
    
    Parameters
    ----------
    dbpath : str
        Database file path.

    Optional Parameters
    -------------------
    start : str
        Pandas-parseable date-time string; default will query from start of database if None.
    end : str
        Pandas-parseable date-time string; default will query from start of database if None.
    tags : list
        List of tags to import; default is None and imports all tags.
    condition : list of strings
        Additional SQL language to specify certain constraints, given as
        a list of conditionals.
    addquery : str
        Additional query language to specify conditions on what data to
        return; default=None.
    tz : str
        The time zone used for the start/end values.
    verbose : bool
        Report additional information during execution.

    Returns
    -------
    Pandas dataframe with query results.
    """
    # convert datetimes to unix for query, or default to limits
    if start is None:
        start = '(SELECT MIN("date_time") from "{}")'.format(table)
    else:
        start = int(pd.Timestamp(start, tz=tz).value/1e9)
    if end is None:
        end = '(SELECT MAX("date_time") from "{}")'.format(table)
    else:
        end = int(pd.Timestamp(end, tz=tz).value/1e9)

    if tags is None:
        cols = '*'
    else:
        cols = '"date_time",'
        for t in tags:
            cols += '"{}",'.format(t)
        cols = cols[:-1]  # delete trailing comma

    query = 'SELECT {cols} FROM "{table}" WHERE "date_time" BETWEEN {start} AND {end}'\
            .format(cols=cols, table=table, start=start, end=end)
    if condition is not None:
        warn("'condition' will be replaced with 'addquery' in the future",
             category=DeprecationWarning)
        for c in condition:
            query += ' AND '+c
    if addquery is not None:
        query += addquery
    
    if verbose:
        print(query)
    
    timeit = Timer()
    engine = create_engine('sqlite:///'+dbpath)
    data = pd.read_sql(query, engine, index_col='date_time', parse_dates=['date_time'])
    # convert back to requested time zone and strip tz-awareness
    data.index = data.index.tz_localize('utc') \
                .tz_convert(tz) \
                .tz_localize(None)
    
    timeit.split('data loaded', report=verbose)
    
    return data


def sqlite_structure_query(dbpath, print_results=True):
    """
    Determine the table and column structure for a db.    
    Parameters
    ----------
    dbpath : str
        Database file path.
    Returns
    -------
    Dictionary of tables and columns.
    """
    conn = sqlite3.connect(dbpath)
    c = conn.cursor()
    c.execute('SELECT name FROM sqlite_master WHERE type="table"')
    tables = set([t[0] for t in c.fetchall()])
    tabcols = {}
    for tab in tables:
        print(tab)
        c.execute('PRAGMA table_info("{}")'.format(tab))
        tabcols[tab] = [[i[1], i[2]] for i in c.fetchall()]
        if print_results:
            for col in tabcols[tab]:
                print('    ', col)
    return tabcols

