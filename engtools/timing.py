#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for timing code execution
"""


import time

class Timer(object):
    """
    A timing object for code execution evaluation.
    """

    def __init__(self):
        """
        Initialize the object with the current time.
        """

        self.t = [time.time()]
        self.lastsplit = 0.0
        self.lasttotal = 0.0


    def _tick(self):
        self.t.append(time.time())
        self.lastsplit = self.t[-1]-self.t[-2]
        self.lasttotal = self.t[-1]-self.t[0]


    def split(self, text='', report=True):
        """
        Record the current time and report the `split` time.

        Optional Parameters
        -------------------
        text : str
            Message to display with the split time.
        report : bool
            Report time if no text was supplied.
        """

        self._tick()
        if report:
            if text != '':
                text += ' '  # add space to maintain format
            print('-------- {:.1f} s {}--------'.format(self.lastsplit, text))


    def total(self, text=''):
        """
        Record the current time and report the overall elapsed time.
        """

        self._tick()
        if text != '':
            text += ' '
        print('-------- total: {:.1f} s {}--------'.format(self.lasttotal, text))
