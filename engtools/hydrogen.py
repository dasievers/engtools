#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hydrogen specific items
"""


import numpy as np
from pandas import read_csv
import os


ppath = os.path.split(__file__)[0]

# =============================================================================
#
# =============================================================================
# check for coefficient table, then import
fname = 'h2_dens_coeffs.csv'
fpath = os.path.join(ppath, fname)
if not os.path.isfile(fpath):
    raise IOError(f'cannot find {fname}')
coeffs = read_csv(fpath)
coeffs.set_index('i', drop=True, inplace=True)


class Hydrogen:
    """
    A hydrogen state. Initialized by temperature and pressure.

    ====== ==================================== =======
    param  description                          unit
    ====== ==================================== =======
    P      pressure                             (MPa)
    T      temperature                          (C)
    ====== ==================================== =======

    Attributes
    ----------
        P : array or value
            The absolute pressure of the hydrogen in (MPa).

        T : array or value
            The absolute temperature of the hydrogen in (K).

        Z : array or value
            The compressibility factor.

    """

    def __init__(self, P, T):
        # initialize and set state
        self._P = P
        self._T = T

        self.calc()

    def __str__(self):
        return f'{self.Z:.3f}'

    def calc(self):
        coefcalc = lambda a, b, c: a * (100/self.T)**b * self.P**c
        items = []
        for i in range(1,10):
            items.append(coefcalc(coeffs.loc[i, 'ai'],
                                  coeffs.loc[i, 'bi'],
                                  coeffs.loc[i, 'ci']))

        items = np.array(items)
        self.Z = 1 + np.sum(items)


    # update dynamic properties
    @property
    def T(self):
        return self._T
    @T.setter
    def T(self, T):
        self._T = T
        self.calc()

    @property
    def P(self):
        return self._P
    @P.setter
    def P(self, P):
        self._P = P
        self.calc()


# =============================================================================
# test
# =============================================================================
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    h = Hydrogen(0.1, 300)

    p_psi = np.linspace(0, 16000)
    p_mpa = p_psi / 14.696 * 101.325 / 1000
    z = []
    for p in p_mpa:
        h.P = p
        z.append(h.Z)
    z = np.array(z)

    fig, ax = plt.subplots()
    ax.plot(p_psi, z)
    ax.set_xlabel('pressure (psi)')
    ax.set_ylabel('Z')




