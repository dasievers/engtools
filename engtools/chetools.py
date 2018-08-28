#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from numpy import exp
from pandas import read_csv
from scipy.interpolate import interp1d

ppath = os.path.split(__file__)[0]

Patm_NREL = 24.02/29.92*101.325  # kPa, atmospheric, average of NREL's SRRL weather station

class SatSteam:
    """
    A saturated steam state. Initialized by a corresponding pressure,
    temperature, etc. State is internally stored as an absolute
    pressure in kPa. Default Patm reference is at NREL's STM site.
    
    Valid parameters are as follows
    
    ====== ============================================
    param   unit
    ====== ============================================
    P       (kPa)
    Pg      (kPa gauge)
    psig    (psi gauge), for pilot-plant compatability
    T       (C)
    v       (m3/kg)
    r       (kg/m3)
    hf      (kJ/kg)
    hg      (kJ/kg)
    hfg     (kJ/kg)
    s       (kJ/kg-K)
    ====== ============================================
    
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
    fpath = os.path.join(ppath, 'steam_table_sat.txt')
    if not os.path.isfile(fpath):
        raise IOError('cannot find steam_table_sat.txt')
    steam = read_csv(fpath, sep='\t', skiprows=(0,1,3))
    
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
        


def henry_constant(T, gas):
    """
    Calculate henry constant for a particular gas at a given temperature.
    
    Parameters
    ----------
    T : float, (C)
        Temperature
    gas : str
        Name of gas.
        
        ======================= =
        available gases
        ======================= =
        oxygen
        ======================= =

    Returns
    -------
    Henry constant in (mM/atm)
    """
    # data from http://www.mpch-mainz.mpg.de/~sander/res/henry.html
    #TODO add all the gases out of the reference by using an external data table
    T = float(T + 273)
    kHstd = {'oxygen':1.3}  # mM/atm
    ddt = {'oxygen':1500}  # K
    return kHstd[gas] * exp(ddt[gas] * (1/T - 1/298.))  # mM/atm


      
    
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









