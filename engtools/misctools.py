#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import re
from datetime import datetime


# =============================================================================
# Helper Functions
# =============================================================================
def smooth(x, window_len=10):
    """
    Adapted from: http://scipy.org/Cookbook/SignalSmooth
    Uses 'flat' window size to return moving average
    
    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    w = np.ones(window_len,'d')
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]


# =============================================================================
# Data Filters
# =============================================================================
def _shortdate(s, df):
    """
    Filter datetimes in shorthand (time only) and add a date to them.
    Useful when a df only occurs over one date and you just need to
    filter over times.
    Uses the very first date of the passed dataframe as the "date" for all, 
    if no date has been found in the passed string.
    Return a pd.datetime object either way.
    """
    hasdate = any([
                    re.match(r'[0-9]{4}', s),  # look for a 4 digit year
                    re.match(r'[0-9]+-[0-9]+', s),  # look for numbers separated by dash
                    re.match(r'[0-9]+/[0-9]+', s),  # look for numbers separated by slash
                    ])
    if hasdate:
        return pd.to_datetime(s)
    else:
        assert pd.to_datetime(df.index[0]).date() == pd.to_datetime(df.index[-1]).date(),\
                'dataframe spans multiple dates'
        d = df.index[0]
        return pd.to_datetime(datetime.combine(pd.to_datetime(d).date(), 
                                               pd.to_datetime(s).time()))
    

def dtrange(df, start=None, end=None):
    """
    Return a dataframe between the two given date/times, even if there is no
    exact match within the df index.
    
    Parameters
    ----------
    df : pandas DataFrame
        Input.
    
    Optional Parameters
    -------------------
    start : str
        Date-time in a format that can be automatically parsed by pandas.
    end : str
        Date-time in a format that can be automatically parsed by pandas.    


    Returns
    -------
    pandas DataFrame between start and end indices
    stard index
    end index
    """
    index = pd.to_datetime(df.index)
    if start is not None:
        startt = _shortdate(start, df)
        idx_start = index.searchsorted(startt)
    else:
        idx_start = 0
        startt = index[0]
    if end is not None:
        endt = _shortdate(end, df)
        idx_end = index.searchsorted(endt)
        df = df.iloc[idx_start:idx_end, :]
    else:
        df = df.iloc[idx_start:, :]
        endt = index[-1]
    return df, startt, endt



def df_smooth(df, win):
    """
    Smooth each series within a pd.DataFrame. If gaps in time larger than
    the window are found, the dataframe is split and separately smoothed
    before recombining.
    
    Index type of datetime64[ns] is automatically detected and the window units
    are in seconds; otherwise window and df index units must match.
    
    Parameters
    ----------
    df : pandas DataFrame
        Input.
    win : float
        Window to smooth over, s for datetime64[ns].
        
    Returns
    -------
    pandas DataFrame containing the converted data
    
    """
    diff = np.diff(df.index.values.astype(float))
    if str(df.index.dtype)=='datetime64[ns]':
        diff = diff/1e9  # s
    # calcualte average dt as the median (avoid large gaps influencing result)
    avgdt = np.median(diff)
    window = int(win//avgdt)  # index window
    
    # test df for time gaps and split up if gaps are greater than the window
    # this is to avoid smoothing over the gaps
    gapidx = np.nonzero(diff > win)[0] + 1
    gapidx = np.insert(gapidx, 0, 0)
    gapidx = np.append(gapidx, df.shape[0])
    print('{} contiguous datasets detected: now smoothing...'.format(len(gapidx)-1))
    splits = []
    for i in range(len(gapidx)-1):
        i1 = gapidx[i]
        i2 = gapidx[i+1]
        splits.append(df.iloc[i1:i2, :])
    # individually smooth each sub-df
    smoothfunc = lambda x: smooth(x, window)
    print('|'+(len(gapidx)-3)*'-'+'|')
    for j in splits:
        print('*', end='')
        if j.shape[0] > window:  # only smooth if slice is large enough
            #TODO need to speed this step up!!
            mask = j.dtypes==float  # only smooth floats
            j.loc[:,mask] = j.loc[:,mask].apply(smoothfunc)
            
            ## old
            # for i in range(len(j.dtypes)):
            #     if 'float' in str(j.dtypes.values[i]).lower():
            #         # only smooth floats
            #         idx = j.dtypes.index[i]
            #         j[idx] = smooth(j[idx].values, window)
                    
    df = pd.concat(splits)
    return df


def indexconvert(df, units, start=0, chopgaps=False, gapthresh=None, dropold=True):
    """
    Change the index of a dataframe to a different unit.  Index type of 
    datetime64[ns] is assumed.
    
    Parameters
    ----------
    df : pandas DataFrame
        Input.
    units : str
        Desired index unit:
            's'
            'min'
            'h'
            'd'
    
    Optional Parameters
    -------------------
    start : float
        t-zero value for the dataset, added to the index.
    chopgaps : bool
        Remove missing time gaps from dataset; requires gapthresh.
    gapthresh : float
        The threshold for removal of time gaps, units same as units.
    dropold : bool
        Drop old index or keep as a new col.

    Returns
    -------
    pandas DataFrame containing the converted data
    
    """
    assert str(df.index.dtype)=='datetime64[ns]'
    df.sort_index(inplace=True)  # ensure sorted by datetime
    
    if units == 's':
        c = 1
    elif units == 'min':
        c = 1./60
    elif units == 'h':
        c = 1./3600
    elif units == 'd':
        c = 1./3600/24
    
    if not chopgaps:
        # just convert the index
        times = c * df.index.values.astype(float)/1e9
        df.reset_index(drop=dropold, inplace=True)
        df.index = times - times[0]
        return df
    
    if chopgaps:    
        diff = c * np.diff(df.index.values.astype(float))/1e9
        # test df for time gaps and split up if gaps are greater than the threshold
        gapidx = np.nonzero(diff > gapthresh)[0] + 1
        gapidx = np.insert(gapidx, 0, 0)
        gapidx = np.append(gapidx, df.shape[0])
    
        splits = []
        for i in range(len(gapidx)-1):
            i1 = gapidx[i]
            i2 = gapidx[i+1]
            splits.append(df.iloc[i1:i2, :])
        # individually convert index for each sub-df
        t0 = start
        for j in splits:
            times = c * j.index.values.astype(float)/1e9  # convert
            times = times - times[0] + t0  # reset reference
            try:
                t0 = times[-1] + (times[-1] - times[-2])  # slightly displaced to avoid duplicate indices
            except IndexError:  # too short of slice
                t0 = times[-1] + 0.0001
            j.reset_index(drop=dropold, inplace=True)
            j.index = times
        dfout = pd.concat(splits)
        return dfout


#==============================================================================
# Local Minima/Maxima Finders
#==============================================================================

def local_minima(x, mode='simple', r=100, c=None):
    '''
    Find local minima in an array and return indices.
    
    Parameters
    ----------
    x : 1D array    

    Optional Parameters
    -------------------
    mode : Operation mode:
            'simple' finds minima via derivative method.
            'width' searches within that +/- range for the minimum at each simple seed.
            'ceiling' defines the local valley by the left and right ceiling and finds minima.
    r : int, +/- indices range to look back and forward for local valley
        (required for 'width' mode)
    c : float, cutoff ceiling value; definition of local valley
        (required for 'ceiling' mode)

    Returns
    -------
    Array with local minima indices of original array, sorted.    
    '''
    # calculate initial guesses by simple differentiation, using gradient to return equal size array 
    g = np.nonzero(np.concatenate((np.diff(np.sign(np.gradient(x)))>0, [False])))[0]
    result = []
    if mode == 'width':
        # for each guess, find minimum value within the +/- lookaround range (r)
        for i in range(len(g)):
            start = g[i]-r
            if start < 0:
                start = 0
            end = g[i]+r
            if end > len(x)-1:
                end = len(x)-1
            result.append(np.argmin(x[start:end]) + start)
    if mode == 'ceiling':
        # for each guess, find minima for local valley under the ceiling cutoff (c)
        # delete g where x[g] > c
        g = np.delete(g, np.nonzero(x[g]>c)[0])
        # go backwards and forwards for each g until x > c and report indices
        for i in range(len(g)):
            idx = g[i]
            while x[idx] < c and idx > 0:
                idx -= 1
            start = idx
            idx = g[i]
            while x[idx] < c and idx < len(x)-1:
                idx += 1
            end = idx
            result.append(np.argmin(x[start:end]) + start)
    result = list(set(result))  # delete duplicates
    result = np.sort(np.array(result))  # sort
    if mode == 'simple':
        # just return the initial derivative values
        result = g
    return result
    

def local_maxima(x, mode='simple', r=100, c=None):
    '''
    Find local maxima in an array and return indices.
    
    Parameters
    ----------
    x : 1D array    

    Optional Parameters
    -------------------
    mode : Operation mode:
            'simple' finds maxima via derivative method.
            'width' searches within that +/- range for the maximum at each simple seed.
            'floor' defines the local valley by the left and right ceiling and finds minima.
    r : int, +/- indices range to look back and forward for local peak
        (required for 'width' mode)
    c : float, cutoff floor value; definition of local peak
        (required for 'floor' mode)

    Returns
    -------
    Array with local minima indices of original array, sorted.    
    '''
    # calculate initial guesses by simple differentiation, using gradient to return equal size array
    g = np.nonzero(np.concatenate((np.diff(np.sign(np.gradient(x)))<0, [False])))[0]
    result = []
    if mode == 'width':
        # for each guess, find maximum value within the +/- lookaround range (r)
        for i in range(len(g)):
            start = g[i]-r
            if start < 0:
                start = 0
            end = g[i]+r
            if end > len(x)-1:
                end = len(x)-1
            result.append(np.argmin(x[start:end]) + start)
    if mode == 'floor':
        # for each guess, find maxima for local peak above the floor cutoff (c)
        # delete g where x[g] > c
        g = np.delete(g, np.nonzero(x[g]<c)[0])
        # go backwards and forwards for each g until x > c and report indices
        for i in range(len(g)):
            idx = g[i]
            while x[idx] > c and idx > 0:
                idx -= 1
            start = idx
            idx = g[i]
            while x[idx] > c and idx < len(x)-1:
                idx += 1
            end = idx
            result.append(np.argmax(x[start:end]) + start)
    result = list(set(result))  # delete duplicates
    result = np.sort(np.array(result))  # sort
    if mode == 'simple':
        # just return the initial derivative values
        result = g
    return result
    
    
#==============================================================================
# Polynomial Fitting Functions -  for surface plots and modeling
#==============================================================================

#def polyfit2d(x, y, z, order=2):
#    # http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
#    import itertools
#    x = np.asarray(x)
#    y = np.asarray(y)
#    z = np.asarray(z)
#    ncols = (order + 1)**2
#    G = np.zeros((x.size, ncols))
#    ij = itertools.product(range(order+1), repeat=2)
#    for k, (i,j) in enumerate(ij):
#        G[:,k] = x**i * y**j
#    m = np.linalg.lstsq(G, z)[0]
#    return m.reshape(n+1,n+1)


def polyfit2d(x, y, f, deg, mask=None):
    '''
    Given coordinates and values, fit a polynomial and return parameters.
    
    Parameters
    ----------
    x : 1D array of coordinates
    y : 1D array of coordinates
    f : 1D array of expression results / data
    deg : degree of polynomial to fit

    Optional Parameters
    -------------------
    mask : 2D array of booleans to exclude factors of full polynomial

    Returns
    -------
    2D array with polynomial factors.
    '''
    from numpy.polynomial import polynomial
    import numpy as np
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    try:
        len(deg)
    except:
        deg = [deg,deg]
    deg = np.asarray(deg)
    vander = polynomial.polyvander2d(x, y, deg)
    # apply mask to delete terms if desired
    if mask != None:
        mask = np.asarray(mask).flatten()
        vander = vander.transpose()[mask].transpose()
#    vander = np.delete(vander, (7,8,9), 1)
#    vander = vander.reshape((-1,vander.shape[-1]))
#    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f)[0]
    # insert zeroes to deleted terms if a mask was applied to keep matrix order
    if mask != None:
        cc = np.zeros(mask.size)
        j = 0
        for i, item in enumerate(mask):
            if item:
                cc[i] = c[j]
                j += 1
    else:
        cc = c
    return cc.reshape(deg+1)
    
    
def polyfit3d(x, y, z, f, deg, mask=None):
    '''
    Given coordinates and values, fit a polynomial and return parameters.
    
    Parameters
    ----------
    x : 1D array of coordinates
    y : 1D array of coordinates
    z : 1D array of coordinates
    f : 1D array of expression results / data
    deg : degree of polynomial to fit

    Optional Parameters
    -------------------
    mask : 3D array of booleans to exclude factors of full polynomial

    Returns
    -------
    3D array with polynomial factors.
    '''
    from numpy.polynomial import polynomial
    import numpy as np
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    f = np.asarray(f)
    try:
        len(deg)
    except:
        deg = [deg,deg,deg]
    deg = np.asarray(deg)
    vander = polynomial.polyvander3d(x, y, z, deg)
    # apply mask to delete terms if desired
    if mask != None:
        mask = np.asarray(mask).flatten()
        print('old', vander.shape)
        vander = vander.transpose()[mask].transpose()
        print('new', vander.shape)
#    vander = vander.reshape((-1,vander.shape[-1]))
#    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f)[0]
    # insert zeroes to deleted terms if a mask was applied to keep matrix order
    if mask != None:
        cc = np.zeros(mask.size)
        j = 0
        for i, item in enumerate(mask):
            if item:
                cc[i] = c[j]
                j += 1
    else:
        cc = c
    return cc.reshape(deg+1)




    
    