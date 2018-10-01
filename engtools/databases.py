#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Database management tools
"""

import pymssql
import pandas as pd
import os
import glob
import numpy as np
import sqlite3
from sqlalchemy import create_engine
from datareading import batch_read
try:
    from timing import Timer
except ModuleNotFoundError:
    # Timer is helpful, but not necessary; create dummy class if not exist
    class Timer:
        def dummy(*args, **kwargs): pass
        def __getattr__(self, _): return self.dummy
from warnings import warn


ppath = os.path.split(__file__)[0]

# =============================================================================
# Ignition Microsoft SQL    
# =============================================================================

def ignition_query(start_time, end_time, taglist, login, addquery=None, verbose=False):
    """
    Query the IBRF Pilot Plant SCADA (Ignition) SQL server 
    for the desired time range and tags. Note, query times are local
    and resulting date-times are local (America/Denver).
    
    Parameters
    ----------
    start_time : str
        Date-time string, parseable by pandas.
    end_time : str
        Date-time string, parseable by pandas.
    taglist : list
        A list of all desired tags. Must be in same format as
        found in SQL tables created by Ignition (consult your own documentation).
    login : list or str
        Server login info: [server, database, user, passwd]. Script will
        automatically find this in user-provided file inside the module
        directory if a string with the name of this file is passed instead.
        The file format is the same items in the above list in ASCII
        document with line break between each item.

    Optional Parameters
    -------------------
    addquery : str
        Additional query language to specify conditions on what data to
        return; default=None.
    verbose : bool
        If True, prints progress & statistics during execution; default=False.

    Returns
    -------
    Pandas DataFrame containing the requested data.
    """

    timeit = Timer()

    # =============================================================================
    # Setup
    # =============================================================================
    # SQL server login
    if type(login) is str:
        # find local login info
        with open(os.path.join(ppath, login), 'r') as f:
            logintext = f.read()
        server, database, user, password = logintext.split('\n')
    
    # join tag paths copied from Ignition tag browser
    assert type(taglist) is list
    tags = r"('" + r"', '".join(taglist) + "')"
    
    # SQL times in UNIX format (ms)
    # must convert given local time to UTC, since database timestamps are in UTC
    epoch_start_time = pd.Timestamp(start_time, tz='America/Denver').value/1e9*1000
    epoch_end_time = pd.Timestamp(end_time, tz='America/Denver').value/1e9*1000
    
    # =============================================================================
    # Query to determine what database partitions to look in for specified range
    # =============================================================================
    # select partitions at the start, end, and everything in-between
    query1 = "SELECT pname FROM sqlth_partitions WHERE \
            ({start} BETWEEN start_time AND end_time) \
            OR ({end} BETWEEN start_time AND end_time) \
            OR (({start} < start_time) AND ({end} > end_time))"\
            .format(start=str(epoch_start_time), end=str(epoch_end_time))
    if verbose:
        print(query1)
    # http://pymssql.org/en/stable/pymssql_examples.html#using-the-with-statement-context-managers
    with pymssql.connect(server, user, password, database) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query1)
            partitions = pd.DataFrame(cursor.fetchall(), columns=['pname'])
    partitions = partitions['pname'].tolist()
    
    if verbose:
        print(partitions)
    timeit.split('database partitions queried for locations', report=verbose)
    # exit early if empty result
    if len(partitions) == 0:
        print('*** ignition_query alert: No data exists for the given date-time range ***')
        return
    
    # =============================================================================
    # Query the partition tables for the tags and range of interest
    # =============================================================================
    query2 = ""
    for i, partition in enumerate(partitions):
        query2 += "SELECT t_stamp, tagpath, floatvalue, intvalue \
                   FROM {partition} JOIN sqlth_te ON tagid = id \
                   WHERE tagpath IN {tag} AND (t_stamp >= {start}) AND (t_stamp <= {end})"\
                  .format(partition=partition, tag=tags, 
                          start=str(epoch_start_time), end=str(epoch_end_time))
        
        if addquery is not None:
            query2 += addquery 
        if i<(len(partitions)-1):
            query2 += " UNION "
        else:  # finish query
            query2 += " ORDER BY t_stamp"
    if verbose:
        print(query2)
    with pymssql.connect(server, user, password, database) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query2)
            data = pd.DataFrame(cursor.fetchall(), columns=['t_stamp', 'tagpath', 'floatvalue', 'intvalue'])
    
    timeit.split('data queried', report=verbose)

    # =============================================================================
    # Clean up and reformat the returned data
    # =============================================================================
    # exit if no data was returned
    if data.shape[0]*data.shape[1] == 0:
        print('*** ignition_query alert: No results were returned for the given date-time range and tags ***')
        print('Values are only recorded when they change outside their deadbands. Try querying a larger date-time range.')
        return
    
    # Creates a column labeled value and fills it with either the Floatvalue or the Intvalue by Summing together
    datana = data.fillna(0)
    data['value'] = datana['floatvalue'] + datana['intvalue']
    
    #Drops the floatvalue and intvalue columns now that value is the only one needed (may not be necessary)
    data = data.drop(['floatvalue','intvalue'],axis = 1)

    # Note that if querying the same data from ignition the times may be 1 second off
    # since it appears ignition converts to an int instead of rounding to the nearest second.
    # We will be rounding to the nearest second for more accurate data.
    # Time is also converted back from ms to s.
    data['t_stamp'] = pd.to_datetime(round(data['t_stamp'].astype(float)/1000), unit='s', utc=True)
    
    data.set_index('t_stamp', drop=True, inplace=True)
    data.index = data.index.tz_convert('America/Denver')
    data.index = data.index.tz_localize(None)
    data['value'] = data['value'].astype(float)
    
    # reconfigures dataframe to have a timestamp index with tagpath columns fills 
    # in the values with floatvalue which now includes the int values.
    # Data assigned to the nearest second (rouded earlier)
    data = data.pivot_table(index=data.index, columns='tagpath', values='value')

    # first forwardfill in NaN values with last, then backfill the beginning of the data if needed
    data.fillna(method='ffill', inplace=True)
    data.fillna(method='bfill', inplace=True)
    
    # check for missing tags (nothing recorded in the queried time span)
    try:
        requested = set(taglist)
        result = set(data.columns)
        missing = requested - result
        if len(missing) > 0:
            print('WARNING: the following tags are missing from the SQL import')
            for i in missing:
                print(i)
            print('This is possibly caused by the date range too narrow to catch')
            print('a recorded value---occurring when a value is not changing.')
    except KeyError:
        pass

    timeit.split('post-processing complete', report=verbose)
    
    return data


class IgnitionTags:
    """
    Query the Ignition SQL server 
    for the current list of tags.
    
    Attributes
    ----------
    taginfo : pd.DataFrame
        A table containing the tag list, including creation and retirement
        dates and tagpaths. Tagpath is also parsed to create an Opto22-like
        tag column, although many of these either don't exist, or are not
        valid names in Opto22.
    """
    def __init__(self, login):
        # SQL server login
        if type(login) is str:
            # find local login info
            with open(os.path.join(ppath, login), 'r') as f:
                logintext = f.read()
            server, database, user, password = logintext.split('\n')

        query = "SELECT tagpath,created,retired FROM sqlth_te"
        with pymssql.connect(server, user, password, database) as conn:
            with conn.cursor() as c:
                c.execute(query)
                result = c.fetchall()
                
        taginfo = pd.DataFrame(result, columns=['tagpath','created','retired'])
        taginfo['created'] = pd.to_datetime(taginfo['created'], unit='ms')
        taginfo['retired'] = pd.to_datetime(taginfo['retired'], unit='ms')
        taginfo['tag'] = taginfo['tagpath'].apply(lambda x: os.path.split(x)[1].upper())
        self.taginfo = taginfo
        

# =============================================================================
# Local Databases (sqlite, mysql, etc)
# =============================================================================
def dbupdate_sqlite(db, datadir, tz='US/Mountain', dropcols=None,
                    selectcols=None, depnames={}):
    """
    Update a sqlite database file with new data or create new if none exists.
        
    Parameters
    ----------
    db : str
        Database file path. 
    datadir : str
        Path of the directory housing the files to be imported. Required if 
        db is a path. Caution, all files will be scanned so only put compatible
        files in this directory.

    Optional Parameters
    -------------------
    tz : str
        The time zone used for the date/time values in data files.
    dropcols : list
        A list of columns to ignore in the source data files.
    selctcols : list
        A list of columns to use; default is all, besides dropcols.
    depnames : dict
        A dictionary used for updating/renaming the column names prior
        to database insert.
    """
    datapath_list = np.asarray(glob.glob(os.path.join(datadir, '*')) )
    # check for new data
    # data will only be imported if not in the database already
    # 'import_history' table contains imported files and their date/time ranges
    
    # import previous file import history and existing column names (if exists)
    if os.path.isfile(db):
        history = set()
        with sqlite3.connect(db) as conn:
            c = conn.cursor()
            c.execute('SELECT "import_history"."file" FROM "import_history"')
            history = set(np.asarray(c.fetchall()).flatten())
            c.execute('PRAGMA table_info("values")')
            exist_cols = {i[1] for i in c.fetchall()}
        # delete entries from datapath_list (globbed) if also exist in history
        # print(datapath_list)
        newmask = np.ones(datapath_list.shape, dtype=bool)
        for i, f in enumerate(datapath_list):
            if os.path.split(f)[1] in history:
                newmask[i] = False
        datapath_list = datapath_list[newmask]
        # print(datapath_list)
    else:
        exist_cols = set()
    # Import New Data
    if len(datapath_list) == 0:
        raise Exception('*** no new files were found to import ***')
    timeit = Timer()    
    # import the data
    masterdata = batch_read(datapath_list, dropcols=dropcols, 
                selectcols=selectcols, depricated_names=depnames, 
                tz=tz)
    
    # make a list of imported files to add to the db history table
    filelist = []
    for file in datapath_list:
        filelist.append(os.path.split(file)[1])
    filedata = pd.DataFrame(filelist, columns=['file'])
        
    timeit.split(report=False)        
    
    # Export to SQL
    # check for new colums and add them if needed
    new_cols = list(set(masterdata.columns.values) - exist_cols)
        
    if os.path.isfile(db):
        with sqlite3.connect(db) as conn:
            c = conn.cursor()
            for i, n in enumerate(new_cols):
                if 'int' in str(masterdata[n].dtype):
                    dtp = 'INT'
                elif 'float' in str(masterdata[n].dtype):
                    dtp = 'FLOAT'
                elif 'str' in str(masterdata[n].dtype):
                    dtp = 'VARCHAR(255)'
                else:
                    dtp = ''
                c.execute('ALTER TABLE "values" ADD "{0}" {1}'.format(n, dtp))
    
    # change index to unix time integers (for speed when querying later)
    masterdata.index = (masterdata.index.values.astype(int)/1e9).astype(int)
       
    #TODO add functionality for using MySQL
    #engine = create_engine('mysql+pymysql://root:root@localhost:8889/'+switch)
    engine = create_engine('sqlite:///'+db)
    masterdata.to_sql("values", engine, if_exists='append', index=True, index_label='date_time')
    filedata.to_sql("import_history", engine, if_exists='append', index=False)
    
    timeit.split('elapsed for SQL insert')



def query_sqlite(db, start=None, end=None, tags=None, condition=None, addquery=None, 
                 tz='US/Mountain', verbose=False):
    """
    Query data from a SQLite database.
    
    Parameters
    ----------
    db : str
        Path to the database.

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
        start = '(SELECT MIN("date_time") from "values")'
    else:
        start = int(pd.Timestamp(start, tz=tz).value/1e9)
    if end is None:
        end = '(SELECT MAX("date_time") from "values")'
    else:
        end = int(pd.Timestamp(end, tz=tz).value/1e9)

    if tags is None:
        cols = '*'
    else:
        cols = '"date_time",'
        for t in tags:
            cols += '"{:}",'.format(t)
        cols = cols[:-1]  # delete trailing comma

    query = 'SELECT {cols} FROM "values" WHERE "date_time" BETWEEN {start} AND {end}'\
            .format(cols=cols, start=start, end=end)
    if condition is not None:
        warn("'condition' will be replaced with 'addquery' in the future",
             category=PendingDeprecationWarning)
        for c in condition:
            query += ' AND '+c
    if addquery is not None:
        query += addquery
    
    if verbose:
        print(query)
    
    timeit = Timer()
    engine = create_engine('sqlite:///'+db)
    data = pd.read_sql(query, engine, index_col='date_time', parse_dates=['date_time'])
    # convert back to requested time zone and strip tz-awareness
    data.index = data.index.tz_localize('utc') \
                .tz_convert(tz) \
                .tz_localize(None)
    
    timeit.split('data loaded', report=verbose)
    
    return data


