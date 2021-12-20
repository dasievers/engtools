#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from numpy import exp
from pandas import read_csv
from scipy.interpolate import interp1d

ppath = os.path.split(__file__)[0]

Patm_NREL = 24.02/29.92*101.325  # kPa, atmospheric

# =============================================================================
#
# =============================================================================

def adiabatic_process(x0, x1, xtype, y0, ytype, k=1.4):
    """
    Adiabatic process

    Inputs parameter type (xtype) and initial/final states (x0, x1).
    Outputs final output state (y1) for type (ytype) given initial state (y0).

    Source: ISBN 0-07-238332-1 p328
    """

    if xtype=='T':
        change = x1/x0
    elif xtype=='P':
        change = (x1/x0)**((k-1)/k)
    elif xtype=='v':
        change = (x0/x1)**(k-1)

    if ytype=='T':
        y1 = y0 * change
    elif ytype=='P':
        y1 = y0 * change**(1/((k-1)/k))
    elif ytype=='v':
        y1 = y0 / change**(1/(k-1))

    return y1


class SatSteam:
    """
    A saturated steam state. Initialized by a corresponding pressure,
    temperature, etc. State is internally stored as an absolute
    pressure in kPa. Default Patm reference is Denver.

    Valid parameters are as follows

    ====== =====================================
    param   unit
    ====== =====================================
    P       (kPa)
    Pg      (kPa gauge)
    psig    (psi gauge), for convenience
    T       (C)
    v       (m3/kg)
    r       (kg/m3)
    hf      (kJ/kg)
    hg      (kJ/kg)
    hfg     (kJ/kg)
    s       (kJ/kg-K)
    ====== =====================================

    Attributes
    ----------
        P : array or value
            The absolute pressure of the steam state in (kPa).

    Methods
    -------
        to(param)
            returns the state of the steam in value units specified
            by `param`.
    """

    # check for steam table, then import
    fname = 'steam_table_sat.csv'
    fpath = os.path.join(ppath, fname)
    if not os.path.isfile(fpath):
        raise IOError(f'cannot find {fname}')
    steam = read_csv(fpath, skiprows=(0,1,3))

    def __init__(self, value, parameter, Patm=Patm_NREL):
        # initialize and set state (held as a pressure value, P)
        self.Patm = Patm
        if parameter=='psig':
            value = value/14.696*101.325  # convert to kPa(g)
            parameter='Pg'
        if parameter=='Pg':  # gauge input for convenience
            parameter = 'P'
            value = value + self.Patm

        x2P_func = interp1d(self.steam[parameter], self.steam['P'], kind='cubic')
        self.P = x2P_func(value)

    def __str__(self):
        return '{:.1f} psig'.format(self.to('psig'))

    def to(self, parameter):
        ptype = None
        if parameter=='psig':
            ptype = parameter
            parameter='P'
        if parameter=='Pg':
            ptype = parameter
            parameter='P'
        P2y_func = interp1d(self.steam['P'], self.steam[parameter], kind='cubic')
        if ptype=='psig':
            return (self.P - self.Patm)/101.325*14.696
        if ptype=='Pg':
            return self.P - self.Patm
        else:
            return P2y_func(self.P)


class Water:
    """
    Water density and viscosity lookup object. Set temperature and read
    results as attributes.

    Attributes
    ----------
        T : array or value
            Temperature, C.

        rho : array or value
            Density, g/cm3.

        vd : array or value
            Dynamic viscosity, N-s/m2.

        vk : array or value
            Dynamic viscosity, m2/s.
    """

    ### density conversions
    # check for data table, then import
    fname = 'density_water_table.csv'
    fpath = os.path.join(ppath, fname)
    if not os.path.isfile(fpath):
        raise IOError(f'cannot find {fname}')
    dens = read_csv(fpath, skiprows=(0,1,3))

    T2rho_func = interp1d(dens['Temperature'],
                dens['Density'], kind='cubic')


    ### viscosity conversions
    # check for data table, then import
    fname = 'viscosity_water_table.csv'
    fpath = os.path.join(ppath, fname)
    if not os.path.isfile(fpath):
        raise IOError(f'cannot find {fname}')
    visc = read_csv(fpath, skiprows=(0,1,3))
    visc['Kinematic viscosity'] = visc['Kinematic viscosity'] / 1e6

    T2vd_func = interp1d(visc['Temperature'],
                visc['Dynamic viscosity'], kind='cubic')
    T2vk_func = interp1d(visc['Temperature'],
                visc['Kinematic viscosity'], kind='cubic')


    # update dynamic properties
    @property
    def T(self):
        return self._T
    @T.setter
    def T(self, T):
        self._T = T
        self.update()


    def __init__(self, T):
        self.T = T


    def update(self):
        self.rho = self.T2rho_func(self.T)

        self.vd = self.T2vd_func(self.T)
        self.vk = self.T2vk_func(self.T)




def henry_constant(T, gas):
    """
    Calculate henry constant for a particular gas at a given temperature.

    Parameters
    ----------
    T : float, (C)
        Temperature.

    gas : str
        Name of gas.

        =============================
        available gases
        =============================
        oxygen
        nitrogen
        hydrogen
        carbon dioxide
        hydrogen sulfide
        ozone
        ammonia
        methane
        =============================

    Returns
    -------
    Henry constant in (mM/atm)
    """

    # data from http://www.mpch-mainz.mpg.de/~sander/res/henry.html
    fpath = os.path.join(ppath, 'henry_coefficients.csv')
    cfs = read_csv(fpath, skiprows=3)
    cfs.set_index('substance', drop=True, inplace=True)
    # kHstd in M/atm
    # ddt in K

    T = float(T + 273)
    kH = cfs.loc[gas, 'kHstd'] * exp(cfs.loc[gas, 'ddt'] * (1/T - 1/298.))  # M/atm
    return kH * 1e3  # mM/atm




# =============================================================================
# test
# =============================================================================
if __name__ == '__main__':
    s1 = SatSteam(500, 'P', Patm=0)
    s2 = SatSteam(100, 'T')
    s3 = SatSteam(200, 'T')
    print(s1.to('hfg'))
    print(s2)
    print(s3.to('v'))

    w = Water(55)
    print(w.vk)
    print(w.vd)
    print(w.rho)

    print(henry_constant(60, 'hydrogen'))









