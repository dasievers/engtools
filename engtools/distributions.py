#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import math
from scipy.integrate import cumtrapz
from scipy import special
from engtools import indexconvert
from scipy import optimize


#==============================================================================
# gamma distribution
#==============================================================================
def _conversion(x, mean, sigma, skew):
    """
    Common calculations.
    Convert mean, sigma, skew to alpha/beta and offset values normally used
    for the Gamma distribution calculation.
    """
    alpha = 4/skew**2
    beta = 0.5*sigma*skew
    x0 = mean - 2*sigma/skew
    y = (x - x0)/beta
    p = np.zeros(x.size)
    return alpha, beta, y, p


def gammapdf(x, mean, sigma, skew, p0=1):
    """
    Produce the Gamma probability density function given an input array and
    parameters. Inputs are statistical values instead of formal alpha, beta
    normally used for Gamma dist.
    
    Parameters
    ----------
    x : 1D array, or discrete value    
    mean : desired mean of the function
    sigma : desired standard deviation of the function
    skew : desied skew of the function

    Optional Parameters
    -------------------
    p0 : scaling parameter, multiplier of returned array

    Returns
    -------
    Array same size as x.
    
    """
    alpha, beta, y, p = _conversion(x, mean, sigma, skew)
    p[y>0] = p0/beta*np.exp((alpha-1)*np.log(y[y>0]) - y[y>0] - special.gammaln(alpha))
    return p


def gammacdf(x, mean, sigma, skew):
    """
    Produce the Gamma cumulative distribution function given an input array and
    parameters. Inputs are statistical values instead of formal alpha, beta
    normally used for Gamma dist.
    
    Parameters
    ----------
    x : 1D array, or discrete value    
    mean : desired mean of the function
    sigma : desired standard deviation of the function
    skew : desied skew of the function

    Returns
    -------
    Array same size as x.
    
    """
    alpha, beta, y, p = _conversion(x, mean, sigma, skew)
    p[y>0] = special.gammainc(alpha, y[y>0])
    return p


#==============================================================================
# Residence Time Distribution Analysis
#==============================================================================
class RTD:
    """
    Perform residence time distribution calculations given arrays of time and
    signal data.  Signal must be a probability distribution function, not
    cumulative distribution function. 
    Calculates the mean time, variance, skewness, and kurtosis
    of the dataset.  Dataset is first normalized so area = 1.
    
    Parameters
    ----------
    OPTION 1: pass separate time array with separate signal array
        t : 1D array
            Time values.
        s, E, F : 1D array (pick one)
            Data signal values.  Must be same length as t.
    OPTION 2: pass a series
    s_series, E_series, F_series : pandas series (pick one)
        Signal with index as datetime object. Will be digested to t and s.

    Optional Parameters
    -------------------
    mean, sigma, skew : float
        Just pass results to object for further calcs.
        
    time_unit: str
        Specify other time unit other than default 'min'. Only used for
        digesting pandas series inputs.


    Attributes
    ----------
    mean
    var : variance
    sigma : standard deviation
    skew : skewness
    kurt : kurtosis
    s_area : signal normailzation divisor. divide raw signal by this to give normalized E curve.
    E : probability distribution function, s divided by s_area.
    F : cumulative distribution function, same size as t and s input arrays. 
    
    Methods
    -------
    s_process
        Process the "s" signal array. Use kwarg 'baseline' to specify
        Auto (True), give a baseline value (float), or 

    Example
    --------
    >>> rtd = RTD(s_series=timedataframe['response'])
    >>> rtd.s_process()
    >>> rtd.calc_params()
    >>> print(rtd.mean, rtd.sigma, rtd.skew)
    
    """
    def _timeconvert(self, ser, time_unit):
        assert type(ser) == pd.core.series.Series
        # extract series data to arrays
        if type(ser.index) == pd.core.indexes.datetimes.DatetimeIndex:
            return indexconvert(ser, time_unit, start=0)
        else:
            return ser
    
    def __init__(self, **kwargs):
        self.time_unit = 'min'  # default
        self.nstdev = 2  # default, used for uncertainty intervals
        self.s_area = 1.  # default
        # populate with kwarg settings
        allowed_keys = set([
                            't',
                            's',
                            'E',
                            'F',
                            's_series',
                            'E_series',
                            'F_series',
                            'mean',
                            'sigma',
                            'skew',
                            'time_unit',
                            'nstdev',
                            ])
        for k, v in kwargs.items():
            if k in allowed_keys:
                setattr(self, k, v)
            else:
                raise TypeError('{} is an invalid keyword argument'.format(k))

        # immediately extract series information
        if hasattr(self, 's_series'):
            series_time = self._timeconvert(self.s_series, self.time_unit)
            self.t = series_time.index.values
            self.s = series_time.values

        elif hasattr(self, 'E_series'):
            series_time = self._timeconvert(self.E_series, self.time_unit)
            self.t = self.series_time.index.values
            self.E = series_time.values

        elif hasattr(self, 'F_series'):
            series_time = self._timeconvert(self.F_series, self.time_unit)
            self.t = self.series_time.index.values
            self.F = series_time.values

    def _extras(self):
        # calculate extra parameters from basic ones
        self.tbreak = self.mean - 2*self.sigma/self.skew
#        self.mode = (4/self.skew**2 - 1)*0.5*self.sigma*self.skew + self.tbreak
        self.Pe = 8 / (math.sqrt(8 * self.sigma**2 / self.mean**2 + 1) - 1)  # Peclet number
        

    def calcF(self):  # F from E
        self.F = cumtrapz(self.E, self.t, initial=0)
        
    def calcE(self):  # E from F
        self.E = np.gradient(self.F, self.t)
        
    def s_process(self, baseline=False):
        """
        Preprocess signal data "s" into E.
        """
        #subtract baseline (opt)
        if baseline is not False:
            if type(baseline) is str:
                if baseline == 'start':
                    baseline = self.s[0]
                elif baseline == 'end':
                    baseline = self.s[-1]
            self.s = self.s - baseline
        #normalize response
        self.scale = np.trapz(self.s, self.t)  # area under original signal
        self.E = self.s / self.scale
        
    def calc_params(self, **kwargs): 
        """
        Calculate distribution parameters from data.
        Automatically generates E if unavailable.
        """
        if not hasattr(self, 'E'):
            if hasattr(self, 's'):
                self.s_process(**kwargs)
            elif hasattr(self, 'F'):
                self.calcE()
        self.mean = np.trapz(self.t * self.E, self.t)
        #BUG: var and/or skew goes negative and kills downstream calcs if baseline subtract is too big
        self.var = np.trapz((self.t - self.mean)**2 * self.E, self.t)
        self.sigma = np.sqrt(self.var)
        self.skew = 1 / (self.var ** (3/2.)) * np.trapz((self.t - self.mean)**3 * self.E, self.t)
        self.kurt = 1 / (self.var ** (4/2.)) * np.trapz((self.t - self.mean)**4 * self.E, self.t) - 3
        self._extras()
        
    def fit_dist(self, dist='gamma', guess=None, param_bounds=None):
        """
        Fit the E data to a defined distribution "dist".
        param_bounds items must match order of distribution function inputs.
        """
        # first make guess if not provided
        if guess is None:
            if not all([hasattr(self, 'mean'),
                       hasattr(self, 'sigma'),
                       hasattr(self, 'skew'),
                       hasattr(self, 'scale')]):
                self.calc_params()  # auto run
            guess = [self.mean, self.sigma, self.skew, self.scale]
        
        # user defined bounds if available, otherwise set to default
        for j in range(4):
            if param_bounds is None:
                param_bounds = np.full((2,4), np.nan)
            if np.isnan(param_bounds[0,j]):
                if j == 0:
                    param_bounds[0,j] = 0
                elif j == 1:
                    param_bounds[0,j] = 0
                elif j == 3:
                    param_bounds[0,j] = 0
                else:
                    param_bounds[0,j] = -np.inf
            if np.isnan(param_bounds[1,j]):
                param_bounds[1,j] = np.inf
            # check that guesses are within bounds
            if not param_bounds[0,j] < guess[j] < param_bounds[1,j]:
                guess[j] = (param_bounds[1,j] + param_bounds[0,j]) / 2
        
        print("param initial guesses", guess)
        print("param bounds", param_bounds)
        
        # optimize fitment of distribution parameters
        if dist == 'gamma':
            param_bounds = (param_bounds[0,0:4],param_bounds[1,0:4])
            params = optimize.curve_fit(gammapdf, self.t, self.s, p0=guess, bounds=param_bounds)
            perr = np.sqrt(np.diag(params[1])) 
            # final parameters
            self.mean = params[0][0]
            self.sigma = params[0][1]
            self.skew = params[0][2]
            self.scale = params[0][3]
            self.var = self.sigma**2
            # uncertainties ('nstdev' standard deviations)
            self.mean_u = self.nstdev * perr[0]
            self.sigma_u = self.nstdev * perr[1]
            self.skew_u = self.nstdev * perr[2]
            # (error propagation)
            self.var_u = self.var * math.sqrt((2)**2 * (self.sigma_u/self.sigma)**2)
        else:
            raise ValueError("invalid distribution function specified")

        #TODO: warn if results are at the borders of bounds

        self._extras()
        # (error propagation)
        self.tbreak_u = self.tbreak * math.sqrt((1)**2 * (self.mean_u/self.mean)**2 + (1)**2 * (2 * self.sigma_u/self.sigma)**2 + (-1)**2 * (2 * self.skew_u/self.skew)**2)
        self.Pe_u = self.Pe * math.sqrt((1)**2 * (self.mean_u/self.mean)**2 + (-1)**2 * (self.sigma_u/self.sigma)**2)
    
    
    
    
# =============================================================================
# testing
# =============================================================================
if __name__ == '__main__':
    # run frenchpress first
    cond = dtrange(fp.process['column outlet conductivity'], 
                  start='2020-02-28 15:46:06',
                  end='2020-02-28 15:50:21')[0]
    baseline = np.mean(cond.values[:20])
    rtd= RTD(s_series=cond)
    rtd.calc_params(baseline=0.95*baseline)
    rtd.fit_dist(param_bounds=np.array([[0,0,0,0],[5,5,2,100]]))
#    rtd.fit_dist()

    plt.plot(rtd.t, rtd.E)
    plt.plot(rtd.t, gammapdf(rtd.t, rtd.mean, rtd.sigma, rtd.skew))
             

    # 2009 RTD test             
    rtd = RTD(t=tor, s=sor)
    rtd.calc_params()
    print(rtd.mean, rtd.var, rtd.skew, rtd.kurt, rtd.scale)
    das.RTDcalc(tor, sor)[:5]         
             
             