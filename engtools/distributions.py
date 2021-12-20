#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import math
from scipy.integrate import cumtrapz
from scipy import special
from pputils.miscellaneous import indexconvert
from scipy import optimize
import matplotlib.pyplot as plt

pltcolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

#==============================================================================
# gamma distribution
#==============================================================================
def _gammaParamConversion(x, mean, std, skew):
    """
    Convert mean, std, skew to alpha/beta and offset values normally used
    for the Gamma distribution calculation. Distribution array also converted.
    """

    alpha = 4/skew**2
    beta = 0.5*std*skew
    x0 = mean - 2*std/skew
    y = (x - x0)/beta
    p = np.zeros(x.size)
    return alpha, beta, y, p


def gammaPDF(x, mean, std, skew, p0=1):
    """
    Produce the Gamma probability density function given an input array and
    parameters. Inputs are statistical values instead of formal alpha, beta
    normally used for Gamma dist.

    Parameters
    ----------
    x : 1D array, or discrete value (e.g. timeseries to plot distribution over)
    mean : desired mean of the function
    std : desired standard deviation of the function
    skew : desied skew of the function

    Optional Parameters
    -------------------
    p0 : scaling parameter, multiplier of returned array

    Returns
    -------
    Array same size as x.

    """

    alpha, beta, y, p = _gammaParamConversion(x, mean, std, skew)
    p[y>0] = p0/beta*np.exp((alpha-1)*np.log(y[y>0]) - y[y>0] - special.gammaln(alpha))
    return p


def gammaCDF(x, mean, std, skew):
    """
    Produce the Gamma cumulative distribution function given an input array and
    parameters. Inputs are statistical values instead of formal alpha, beta
    normally used for Gamma dist.

    Parameters
    ----------
    x : 1D array, or discrete value (e.g. timeseries to plot distribution over)
    mean : desired mean of the function
    std : desired standard deviation of the function
    skew : desied skew of the function

    Returns
    -------
    Array same size as x.

    """

    alpha, beta, y, p = _gammaParamConversion(x, mean, std, skew)
    p[y>0] = special.gammainc(alpha, y[y>0])
    return p


#==============================================================================
# Residence Time Distribution Analysis
#==============================================================================
class RTD:
    """
    Core residence time distribution calculations given arrays of time and
    signal data.  Signal must be a probability distribution function, not
    cumulative distribution function prior to fitting and calculation.
    Calculates the mean time, variance, skewness, and kurtosis
    of the dataset.  Dataset is first normalized so area = 1.

    Currently, the Gamma distribution is used to model the RTD and a complete
    object will be populated with mean, std, and skew values that completely
    define the RTD. From this state, the object can be used for further analysis
    including comparison or reaction kinetics studies (by other functions).

    Child classes of this class will contain functionality to import and
    preprocess data to enable use of the core functionality here.

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
    name : name of the object for plotting/saving
    mean, std, skew : float (provide all)
        Populate object with results for further calcs.

    time_unit: str
        Specify other time unit other than default 'min'. Only used for
        digesting pandas series inputs.

    Attributes
    ----------
    name : name of the object for plotting/saving
    t : time array
    s : signal array
    mean
    mode
    var : variance
    std : standard deviation
    skew : skewness
    kurt : kurtosis
    scale : signal normailzation divisor. divide raw signal by this to give
            normalized E curve.
    E : probability distribution function, s divided by scale.
    F : cumulative distribution function, same size as t and s input arrays.

    Methods
    -------
    s2E
        Process the "s" signal array. Use kwarg 'baseline' to specify
        Auto (True), give a baseline value (float), or
    calc_params
        Calculate statistical parameters based on raw t & s
    fit_dist
        Fit t & s data to a distribution (gamma) and store statistical
        parameter results.
    save_stats_classic
        Save the statistical parameters to file.
    plt_dist
        Plot the distribution

    Example
    -------
    >>> rtd = RTD(s_series=timedataframe['signal'])
    >>> rtd.s2E()
    >>> rtd.calc_params()
    >>> print(rtd.mean, rtd.std, rtd.skew)
    >>> rtd.fit_dist()
    >>> rtd.save_stats_classic('output.tsv')
    """

    def _timeconvert(self, ser):
        """Convert index, but only if a datetimeindex."""

        assert type(ser) == pd.core.series.Series
        # extract series data to arrays
        if type(ser.index) == pd.core.indexes.datetimes.DatetimeIndex:
            return indexconvert(ser, self.time_unit, start=0)
        else:
            return ser


    # defaults for class
    name = None
    time_unit = 'min'
    nstdev = 2.  # used for uncertainty intervals
    s_area = 1.

    def __init__(self, **kwargs):
        # populate with kwarg settings
        allowed_keys = set([
                            'name',
                            't',
                            's',
                            'E',
                            'F',
                            's_series',
                            'E_series',
                            'F_series',
                            'mean',
                            'std',
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
            series_time = self._timeconvert(self.s_series)
            self.t = series_time.index.values
            self.s = series_time.values

        elif hasattr(self, 'E_series'):
            series_time = self._timeconvert(self.E_series)
            self.t = self.series_time.index.values
            self.E = series_time.values

        elif hasattr(self, 'F_series'):
            series_time = self._timeconvert(self.F_series)
            self.t = self.series_time.index.values
            self.F = series_time.values


    def _extras(self):
        # calculate extra parameters from basic ones
        self.tbreak = self.mean - 2*self.std/self.skew
        self.mode = (4/self.skew**2 - 1)*0.5*self.std*self.skew + self.tbreak
        self.Pe = 8 / (math.sqrt(8 * self.std**2 / self.mean**2 + 1) - 1)  # Peclet number


    #TODO: add decorators for automatic execution when either of these change
    def calcF(self):  # F from E
        self.F = cumtrapz(self.E, self.t, initial=0)


    def calcE(self):  # E from F
        self.E = np.gradient(self.F, self.t)


    def s2E(self, baseline=None):
        """
        Preprocess signal data "s" into E.
        Pass optional 'baseline' as 'start'/'end' point or actual value.
        """

        #subtract baseline (opt)
        if baseline is not None:
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
        """Calculate distribution parameters from data."""

        if not hasattr(self, 'E'):
            self.s2E(**kwargs)
        self.mean = np.trapz(self.t * self.E, self.t)
        #BUG: var and/or skew goes negative and kills downstream calcs if baseline subtract is too big
        # note these eqs are from DAStools and were working in the past (at least thought to)
        self.var = np.trapz((self.t - self.mean)**2 * self.E, self.t)
        self.std = np.sqrt(self.var)
        self.skew = 1 / (self.var ** (3/2.)) * np.trapz((self.t - self.mean)**3 * self.E, self.t)
        self.kurt = 1 / (self.var ** (4/2.)) * np.trapz((self.t - self.mean)**4 * self.E, self.t) - 3
        self._extras()


    def fit_dist(self, dist='gamma', guess=None, param_bounds=None):
        """
        Fit the E data to a defined distribution "dist".
        param_bounds items must match order of distribution function inputs.
        [mean, std, skew, scale]

        Default bounds are 0 - inf, but the skew could be negative
        in certain special situations. In those cases, specify
        bounds explicitly.

        Note about 'meta' attribute:
            Parameter guesses and bounds may be passed as 'fitparams' of the
            meta dictionary attribute. If available, they will be used instead
            of defaults or rough estimations. The 'fitparams' value is another
            dictionary with any of the following keys: mean, std, skew, scale.
            Not all are needed, only the ones that are wish to be set.
            Each of those items has a list with exactly three values:
            [guess, lowerbound, upperbound]. Use np.nan if any of the values
            are not set.
        """

        # set up the guess and param_bounds
        paramlist = ['mean', 'std', 'skew', 'scale']

        # make guess if not explicitly provided
        if guess is None:
            # use initial guesses as calculated parameters from rough data points
            attrpresent = [hasattr(self, p) for p in paramlist]
            if not all(attrpresent):
                self.calc_params()  # auto run
            guess = [getattr(self, p) for p in paramlist]

            # if metadata guesses are also present, override auto guesses
            if hasattr(self, 'meta'):
                if 'fitparams' in self.meta.keys():
                    for i, p in enumerate(paramlist):
                        try:
                            x = self.meta['fitparams'][p][0]  # get first value of list
                            if not np.isnan(x):
                                guess[i] = x
                        except(KeyError):
                            pass  # parameter absent: just skip and leave auto guess value

        # make bounds if not provided
        if param_bounds is None:
            # initialize to defaults
            param_bounds = np.vstack([np.full(len(paramlist), 0),
                                      np.full(len(paramlist), np.inf)])

            # if metadata bounds are also present, override defaults
            if hasattr(self, 'meta'):
                if 'fitparams' in self.meta.keys():
                    for j, p in enumerate(paramlist):
                        try:
                            x = self.meta['fitparams'][p][1]
                            if not np.isnan(x):
                                param_bounds[0,j] = x
                            x = self.meta['fitparams'][p][2]
                            if not np.isnan(x):
                                param_bounds[1,j] = x
                        except(KeyError):
                            pass  # parameter absent: just skip and leave auto guess value

        # check that guesses are within bounds and autocorrect if not
        for j, p in enumerate(paramlist):
            if not param_bounds[0,j] <= guess[j] <= param_bounds[1,j]:
                guess[j] = 1.  # default
                print('Guess is outside parameter bounds for {} and was corrected.'.format(p))

        print("parameters are", paramlist)
        print("param initial guesses", guess)
        print("param bounds", param_bounds)

        # optimize fitment of distribution parameters
        if dist == 'gamma':
            results = optimize.curve_fit(gammaPDF, self.t, self.s,
                                         p0=guess,
                                         bounds=param_bounds)
            perr = np.sqrt(np.diag(results[1]))
            # final parameters
            self.mean = results[0][0]
            self.std = results[0][1]
            self.skew = results[0][2]
            self.scale = results[0][3]
            self.var = self.std**2
            # uncertainties ('nstdev' standard deviations)
            self.mean_u = self.nstdev * perr[0]
            self.std_u = self.nstdev * perr[1]
            self.skew_u = self.nstdev * perr[2]
            # (error propagation)
            self.var_u = self.var * math.sqrt((2)**2 * (self.std_u/self.std)**2)
        else:
            raise ValueError("invalid distribution function specified")
            ### placeholder: add different distributions later on here

        self._extras()
        # (error propagation)
        self.tbreak_u = self.tbreak * math.sqrt((1)**2 * (self.mean_u/self.mean)**2 + (1)**2 * (2 * self.std_u/self.std)**2 + (-1)**2 * (2 * self.skew_u/self.skew)**2)
        self.Pe_u = self.Pe * math.sqrt((1)**2 * (self.mean_u/self.mean)**2 + (-1)**2 * (self.std_u/self.std)**2)

        print(40*'-')
        print(8*'*' + '  FIT REPORT  ' + 8*'*')
        print(self.name)
        print('')
        print(u'breakthrough: {:.2f} {}'.format(self.tbreak, self.time_unit))
        print(u'mean:         {:.2f} \u00B1 {:.2f} {}'.format(
                self.mean, self.mean_u, self.time_unit))
        print(u'std dev:      {:.2f} \u00B1 {:.2f} {}'.format(
                self.std, self.std_u, self.time_unit))
        print(u'skew:         {:.2f} \u00B1 {:.2f} {}'.format(
                self.skew, self.skew_u, self.time_unit))
        print(u'Peclet nr:    {:.1f} \u00B1 {:.1f} {}'.format(
                self.Pe, self.Pe_u, self.time_unit))
        print(40*'-')


    def save_stats(self, fname=None, directory=None):
        """
        Save stats of distribution in text format.

        File name created automatically from object name if
        none is passed.
        Directory is pwd unless specified.
        """

        if fname is None:
            # create from name
            objname = 'RTD'
            if self.name is not None:
                objname = self.name
            fname = objname + ' stats.csv'
        if directory is not None:
            fname = os.path.join(directory, fname)

        stats_key = ['breakthrough', 'mode', 'mean', 'std dev', 'skew', 'Pe']
        stats = [self.tbreak, self.mode, self.mean, self.std, self.skew, self.Pe]
        stats_u = [np.nan, np.nan, self.mean_u, self.std_u, self.skew_u, self.Pe_u]

        df = pd.DataFrame({'stats' : stats,
                           'uncertainty' : stats_u},
                          index=stats_key)
        df.to_csv(fname)


    def plt_dist(self, dist='gamma', plotdots=True, normalized=True, save=False):
        """Plot the fitted distribution E curve with or without data"""

        try:
            tmax = self.meta['runlength']
        except(ValueError, KeyError, AttributeError):
            tmax = max(rtd.t) * 1.15
        tcurve = np.linspace(0, tmax, 500)
        if dist == 'gamma':
            Ecurve = gammaPDF(tcurve, self.mean, self.std, self.skew)
        else:
            raise ValueError("invalid distribution function specified")
        fig, axs = plt.subplots(1,1)
        axs.set_title('RTD for {}'.format(self.name))
        axs.set_xlabel('time ({})'.format(self.time_unit))
        if normalized:
            axs.plot(tcurve, Ecurve, linewidth=2, color=pltcolors[0])
            if plotdots:
                axs.plot(self.t, self.s/self.scale, 'o', markersize=5, color='k')
            axs.set_ylabel('E (1/{})'.format(self.time_unit))
        else:
            axs.plot(tcurve, Ecurve*self.scale, linewidth=3, color=pltcolors[0])
            if plotdots:
                axs.plot(self.t, self.s, 'o', markersize=6, mec='k', mfc='none')
            axs.set_ylabel('signal')

        if save:
            fname = 'RTD curve.pdf'
            if self.name is not None:
                fname = self.name + ' ' + fname
            plt.savefig(fname, bbox_inches='tight')


# =============================================================================
# testing
# =============================================================================
if __name__ == '__main__':
    userdir = os.path.expanduser('~')
    # testing on an RTD run in the FCIC french press apparatus
    from pputils import FrenchPress

    fp = FrenchPress(name='frenchpress_rtd_test')
    fp.load_scada(start='2020-02-28 15:46:06',
                  end='2020-02-28 15:50:21',
                  localdbpath=os.path.join(userdir, 'databases/FCIC/'))
    cond = fp.process['column outlet conductivity']
    baseline = np.mean(cond.values[:40])
    rtd = RTD(s_series=cond)
    rtd.s2E(baseline=0.95*baseline)
    rtd.calc_params()
    rtd.fit_dist()
    rtd.plt_dist(normalized=False)








