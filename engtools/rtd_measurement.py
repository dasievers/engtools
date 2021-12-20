#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RTD measurement subclass, elements
"""


from io import StringIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pputils.miscellaneous import dtrange, indexconvert
from scipy import interpolate as ipt
from scipy import stats

try:
    import regularsmooth as ds  #Jonathn's smoothing algorithm
except ModuleNotFoundError:
    # load module without regsmooth (limited functionality)
    print('NOTICE: regularsmooth could not be imported.\n'\
          +'outlier_filter method will be unavailable.')
    def regularsmooth():
        return None

pltcolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# =============================================================================
#
# =============================================================================
def _simpleaxis(ax):  # removes unnecessary axes from charts
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

# =============================================================================
#
# =============================================================================
from pputils.distributions import RTD
class RTDmeasure(RTD):
    """
    Class for importing and processing RTD data associated with reactors.

    Attributes
    ----------
    *name : name of the object for plotting/saving
    *meta : dictionary of metadata
    *metatable : table of metadata specific to each sub-run of distribution measurement
    signalstream : raw data in pandas Series (after any prefiltering).
        This is user-specified by selecting
        Metadata will later be used to slice these data for actual runs.
    data : pd Dataframe after aligning multiple repeat runs on top
        of each other and converting index to elapsed time (float).
        These are the data used for distribution fitting, etc.

    Methods
    -------
    preprocess_squarewave
        Process the SCADA data for the NaCl/conductivity method that
        needs to be analyzed as a square wave amplitude only valid
        when discharge valves are at certain position (cell filled).
    add_runrawdata
        Add datasets to run dictionary attribute. Can be run multiple
        times to populate with different runs.
    process_runs
        Process raw data signals and populate one Series/Dataframe with
        all the runs. Must either have run preprocess_squarewave or
        add_rundata first.
    outlier_filter
        Remove outlier datapoints either through automatic detection
        or explicit cherrypicking.
    prepare
        Set the t & s attributes from Series or Dataframe objects of data.

    Example
    -------
    >>> rtd = RTDmeasure()
    >>> rtd.load_metadata_classic(metafile='metadata_metsoRTDrunX.txt',
                                  rawdir='/Users/user/RTDstuff')
    >>> rtd.load_process()
    >>> rtd.preprocess_squarewave()
    >>> rtd.process_runs()
    >>> rtd.outlier_filter()
    >>> rtd.prepare()
    >>> rtd.fit_dist()
    >>> rtd.plt_dist()
    >>> rtd.save_stats_classic('output.tsv')
    """

    def _dataloader(self, src, **kwargs):
        if '.xls' in src:
            return pd.read_excel(src, **kwargs)
        else:
            with open(src, 'r', errors='ignore') as f:
                # note errors must be ignored to deal with "can't decode byte 0xff in position 0"
                s = f.read()
            if '\t' in s:
                delimiter = '\t'
            elif ',' in s:
                delimiter = ','
            else:
                delimiter = None  # error
            return pd.read_csv(StringIO(s), delimiter=delimiter, **kwargs)


    def add_rawdata(self, runID, time, signal,
                       file=None, start=None, end=None,
                       **kwargs):
        """
        Add raw data to object. Different functionality depending
        on what is passed:
        1. Add arrays directly:
            Pass runID as str or int
            Pass time as array
            Pass signal as array

        2. Load a file (excel or csv/tsv):
            Pass runID as str or int
            Pass time as str, name of column to use in file
            Pass signal as str, name of column to use in file
            Pass file path

        If time array is datetime (not elapsed time), then
        it is sliced by passed start/end datetimes
        and then index converted to float (unit specified by
        'time_unit' attribute).

        This can be run multiple times to add to the 'rawdata';
        to delete datasets, you must start over with new object.
        Each time this is run, the attribute 'data' is
        made/overwritten as a copy of rawdata.

        Extra kwargs are passed to either pd.read_excel or pd.read_csv,
        depending on the file extension.
        """

        if not hasattr(self, 'rawdata'):
            # initialize
            self.rawdata = pd.DataFrame()

        if file is None:
            # direct import
            s = pd.Series(signal, index=time)
        else:
            # import time/signal from raw data using column names
            df = self._dataloader(file, **kwargs)
            s = pd.Series(df[signal].values, index=df[time].values)
        s.name = 'signal'

        if type(s.index) == pd.core.indexes.datetimes.DatetimeIndex:
            # convert index
            s = dtrange(s, start=start, end=end)
            s = indexconvert(s, self.time_unit)

        s.dropna(how='any', inplace=True)
        s.index.name = f'time ({self.time_unit})'

        rundf = pd.DataFrame(s)
        rundf['run ID'] = runID
        self.rawdata = self.rawdata.append(rundf)
        self.rawdata.sort_index(inplace=True)

        # make copy as 'data' to be worked on
        self.data = self.rawdata.copy(deep=True)
        self.prepare()


    def reset(self):
        """Reset 'data' to 'rawdata'."""

        self.data = self.rawdata.copy()


    def plt_data(self, runID=None, raw=False):
        """
        Plot data or rawdata.
        Specify which run or default all.
        """

        if raw:
            df = self.rawdata
        else:
            df = self.data

        if runID is None:
            title = 'Signal for All Runs'
            s = df['signal']
        else:
            title = f'Signal for Run {runID}'
            s = df[df['run ID'] == runID]['signal']

        fig, axs = plt.subplots()
        axs.plot(s, 'ro')
        axs.set_title(title)
        axs.set_xlabel(f'elapsed time ({self.time_unit})')


    def filter_offset(self, offset, runID=None):
        """
        Add offset to index to accomodate things like delay
        in tracer actually getting to region of interest.
        Pass negative value for delay.
        """

        if runID is None:
            self.data.index += offset
        else:
            idx = self.data['run ID'] == runID
            # make copy, then drop the run
            df = self.data[idx]
            self.data = self.data[~idx]  # keep opposite of what is being worked on

            df.index += offset

            # add back result and re-sort
            self.data = self.data.append(df)
            self.data.sort_index(inplace=True)

        self.prepare()


    def filter_baseline(self, baseline=0, dose=1, runID=None):
        """
        Subtract baseline from signal.
        Then, divide signal by dose.
        Optionally for just one run as specified.
        """

        if runID is None:
            self.data.loc[:, 'signal'] -= baseline
            self.data.loc[:, 'signal'] /= dose
        else:
            idx = self.data['run ID'] == runID
            # make copy, then drop the run
            df = self.data[idx]
            self.data = self.data[~idx]  # keep opposite of what is being worked on

            df.loc[:, 'signal'] -= baseline
            df.loc[:, 'signal'] /= dose

            # add back result and re-sort
            self.data = self.data.append(df)
            self.data.sort_index(inplace=True)


        self.prepare()


    def filter_clip(self, start=None, stop=None, runID=None):
        """
        Clip off start and/or end times of dataset, optionally
        for the specified runID.
        """
        if runID is None:
            self.data = dtrange(self.data, start, stop)
        else:
            idx = self.data['run ID'] == runID
            # make copy, then drop the run
            df = self.data[idx]
            self.data = self.data[~idx]  # keep opposite of what is being worked on

            dfclip = dtrange(df, start, stop)  # clip it

            # add back result and re-sort
            self.data = self.data.append(dfclip)
            self.data.sort_index(inplace=True)

        self.prepare()


    def filter_outliers(self, chauvenet=False, manual=False, plot=True):
        """
        SIGNAL FILTERING / OUTLIER REMOVAL
        Smoothing parameters are established that are used for Chauvenet outlier
        removal and also set up parameters used for constrained smoothing that is
        used for final parameter initial guess (later). Mnaual removal also
        available if indices are passed.
        Duplicate indices from various methods are OK.

        if manual, specify a list of outlier indices to remove (not just True).
        """

        # initialize
        t = self.data.index.values
        s = self.data['signal'].values
        outidx = []

        # Calculate smoothing parameters (used for Chauvenet outlier removal and later).  Code by JJS
        Nhat = 200
        th = np.linspace(0, self.data.index.max(), Nhat)
        try:
            lmbd = self.meta['lmbd_set']
            sh = ds.smooth_data(t, s, xhat=th, lmbd=lmbd)
        except(AttributeError, KeyError):
            sh, lmbd = ds.smooth_data(t, s, xhat=th, lmbd_guess=2.0e-4)
        self.lmbd = float(lmbd) # JJS: bug in regularsmooth if type(lmbd) = np.float64 !!??
        print('JJS smoothing --- lmbd:', self.lmbd)

        # CHAUVENET OUTLIER REMOVAL
        # apply Chauvenet's criteria to remove outliers
        if chauvenet:
            shi = ipt.interp1d(th, sh, 'cubic')(t)
            dels = s - shi # difference between smoothed curve and data
            stdev = np.std(dels,ddof=1)
            # probability allowed for test condition
            prob = 1 - 1./(2*t.size)
            # determine normal distribution value corresponding to 'prob'
            nsd = stats.norm.ppf(1 - (1-prob)/2)
            x = list(np.nonzero(np.abs(dels) > nsd*stdev)[0])
            outidx = outidx + x

        # MANUAL OUTLIER REMOVAL
        # manually pick out bad data using passed indices
        if manual:
            outidx = outidx + manual

        # finish up process
        if plot:
            # check outlier handling
            fig, axs = plt.subplots(1,1)
            axs.set_title('Outlier Analysis')
            axs.plot(th, sh, 'k')
            try:
                axs.plot(th, sh - nsd*stdev, 'r')
                axs.plot(th, sh + nsd*stdev, 'r')
            except(NameError):
                pass
            axs.plot(t[outidx], s[outidx], 'rx', ms=7, mew=2)
            for r in self.data['run ID'].unique():
                axs.plot(self.data['signal'][self.data['run ID'] == r],
                         linestyle='none', marker='$'+str(r)+'$', color='k')
            axs.set_xlabel('time ({})'.format(self.time_unit))
            axs.set_ylabel('signal')

        if len(outidx) == 0:
            print('no outliers were removed')
        else:
            print('outliers removed at indices', outidx)
        # remove outliers from indices
        self.data_dropped = self.data.iloc[outidx, :]
        self.data.drop(self.data.index[outidx], inplace=True)

        self.prepare()


    def prepare(self):
        """
        Set up object for further calculations using RTD base class methods.
        """

        self.t = self.data.index.values
        self.s = self.data['signal'].values






