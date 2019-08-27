#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from engtools import SatSteam
from pputils.plotting import yield_bars
import pandas as pd
import copy

# =============================================================================
# General
# =============================================================================
def fis(ss, ls):
    """
    Calculate no-wash insoluble solids fraction of a slurry.

    Parameters
    ----------
    ss : array or value
        Total solids fraction of whole slurry.
    ls : array or value
        Total solids fraction of slurry liquor.

    Returns
    -------
    Value or array with fis.
    """
    return (ss - ls) / (1 - ls)


def _mol_weights(ctype):
    if ctype=='pentose':
        Mstruct = 132
        Msugar = 150
        Mdeg = 96
    elif ctype=='hexose':
        Mstruct = 162
        Msugar = 180
        Mdeg = 126
    elif ctype is None or ctype=='lignin':
        # nullify MW contributions in calculations; return mass yield instead
        Mstruct, Msugar, Mdeg = 1,1,1
    else:
        raise ValueError('invalid ctype')
    return Mstruct, Msugar, Mdeg


def _molarizer(self, num, den, v, k):
    # for totalizing over component yields
    # adjust for molar basis on the feedstock
    Mstruct, Msugar, Mdeg = _mol_weights(self.ctype[k])
    num += v * self.struct0[k] / Mstruct
    den += self.struct0[k] / Mstruct
    return num, den


class Reactor:
    """
    Base class for reactor type yield calculations.
    """
    def setctype(self):
        """
        Automatically make the ctype dict based on names of feed dicts.
        """
        if hasattr(self, 'struct0'):
            keys = self.struct0.keys()
        elif hasattr(self, 'struct1'):
            keys = self.struct1.keys()
        else:
            raise NameError('structural dictionary not defined; can\'t set ctype')
        self.ctype = {}
        for k in keys:
            if 'lig' in k.lower():
                ct = 'lignin'
            elif 'x' in k.lower() or 'ar' in k.lower():
                ct = 'pentose'
            elif 'g' in k.lower() or 'ce' in k.lower() \
            or 'fr' in k.lower() or 'ma' in k.lower():
                ct = 'hexose'
            else:
                raise ValueError('could not automatically set ctype from structural key `{}`'.format(k))
            self.ctype[k] = ct
            
    def calcy(self, func, iterkeys, **kwargs):
        """
        Generic component and overall yield calculator.
        Used for structural carbs and sugars.
        Must return the values and set them on the call.
        """
        yd = {}
        num, den = 0.0, 0.0
        for k in iterkeys:
            # pull out the k'th item for each kwarg value that is a nested dictionary
            kwargk = copy.copy(kwargs)
            for key, value in kwargk.items():
                if type(value) is dict:
                    kwargk[key] = value[k]
            # run the called yield function with the kwargs for k
            y = func(**kwargk)
            # note, `overall yield` does not consider lignin
            yd[k] = y
            if self.ctype[k] != 'lignin':
                num, den = _molarizer(self, num, den, y, k)
        y_overall = num / den
        return yd, y_overall

    def calcsug_abs(self, sugar):
        # calculate absolute yields of sugars in kg component / kg dry feedstock
        yd = {}
        y_overall = 0
        for k, component in sugar.items():
            if self.ctype[k] != 'lignin':
                Mstruct, Msugar, Mdeg = _mol_weights(self.ctype[k])
                y = self.struct0[k]/Mstruct*component*Msugar
                yd[k] = y
                y_overall += y
        return yd, y_overall
    
    def calcmoleclose(self):
        # calculate mole closures by summing up the component yields
        # lignin not included
        self.mc = {}
        for k in self.y_struct.keys():
            if self.ctype[k] != 'lignin':
                if hasattr(self, 'y_deg'):
                    if k in self.y_struct.keys() and k in self.y_total.keys() \
                    and k in self.y_deg.keys():
                        mc = self.y_struct[k] + self.y_total[k] + self.y_deg[k]
                    else:
                        mc = np.nan
                else:
                    if k in self.y_struct.keys() and k in self.y_total.keys():
                        mc = self.y_struct[k] + self.y_total[k]
                    else:
                        mc = np.nan
                self.mc[k] = mc
        if hasattr(self, 'y_deg'):
            self.omc = self.y_struct_overall + self.y_total_overall + self.y_deg_overall
        else:
            self.omc = self.y_struct_overall + self.y_total_overall

    def plot(self, key=None, **kwargs):
        """
        Make a bar plot of the molar yields. Additional kwargs may
        be passed as defined by `yield_bars`.
        
        Parameters
        ----------
        key : str
            The key of the component family. If None, the overall totals
            will be plotted.
        """
        # assign dummies to nonexistent variables
        if hasattr(self, 'y_struct_overall'):
            y_struct_overall = self.y_struct_overall
        else:
            y_struct_overall = np.zeros(len(self.struct0[list(self.struct0.keys())[0]]))
        if hasattr(self, 'y_oligo_overall'):
            y_oligo_overall = self.y_oligo_overall
        else:
            y_oligo_overall = np.zeros(len(self.struct0[list(self.struct0.keys())[0]]))
        if hasattr(self, 'y_mono_overall'):
            y_mono_overall = self.y_mono_overall
        else:
            y_mono_overall = np.zeros(len(self.struct0[list(self.struct0.keys())[0]]))
        if hasattr(self, 'y_deg_overall'):
            y_deg_overall = self.y_deg_overall
        else:
            y_deg_overall = np.zeros(len(self.struct0[list(self.struct0.keys())[0]]))

        if key is None:
            yield_bars(y_struct_overall, y_oligo_overall, 
                       y_mono_overall, y_deg_overall, title='overall',
                       **kwargs)
        else:
            yield_bars(self.y_struct[key], self.y_oligo[key], 
                       self.y_mono[key], self.y_deg[key], title=key, **kwargs)

    def export(self, exportpath):
        """
        Export the entire contents of reactor object to Excel spreadsheet
        for a complete data-dump to share with others.
        Each attribute has its own sheet.
        
        Parameters
        ----------
        exportpath : str
            Path to save the exported data.

        """
        writer = pd.ExcelWriter(exportpath)
        for k, v in self.__dict__.items():
            # print(k)
            if type(v)==dict:
                try:
                    edf = pd.DataFrame(v)
                except ValueError:
                    edf = pd.Series(v, name=k)
                # edf = pd.DataFrame()
                # for kk, vv in v.items():
                #     edf = edf.append(pd.Series(vv, name=kk))
            else:
                try:
                    edf = pd.Series(v, name=k)
                except TypeError:
                    edf = pd.Series([v], name=k)
            # if edf.
            edf.to_excel(writer, sheet_name=k)
        writer.save()




# =============================================================================
# Pretreatment Yield Calculations
# =============================================================================
"""
Notes on the Methods
--------------------
Yields calculated here are on a molar basis and make a key assumption: total solids
flow in the reactor is conserved and measurement of actual inlet and outlet
flow rates are not necessary. They cancel out, and the key flow split
at the flash tank---slurry versus condensed flash vapor---is calculated via
thermodynamics. This is a key `feature`, since inefficiencies of the condensers
and the resulting inaccuracies in condensate flow rate measurement are
inconsequential.
"""

def rx_solids(TSf, Tr, Tf, cps=1.463):
    """
    Estimate the total solids fraction inside the pretreatment reactor just
    prior to flashing, given the outlet stream properties.
    
    Parameters
    ----------
    TSf : array or value
        Total solids fraction slurry after flash.
    Tr : array or value, (C)
        Reactor temperature.
    Tf : array or value, (C)
        Flash tank temperature.
    cps : array or value, (kJ/kg-K)
        Heat capacity of solids (default is corn stover).

    Returns
    -------
    Array same size as struct.
    """
    steam1 = SatSteam(Tr, 'T')
    steam2 = SatSteam(Tf, 'T')
    hf1 = steam1.to('hf')
    hf2 = steam2.to('hf')
    hg2 = steam2.to('hg')
    
    TSr = (hg2 - hf1)/(hf2 - hf1 + (hg2  - hf2)/TSf + cps*(Tr - Tf))
    return TSr
   

def pt_yield_struct(TS, IS, struct, struct0):
    # note, no MW conversions are necessary since MW is equal for substrate and product
    # enforce boundaries
    struct = np.clip(struct, 0, 1)
    struct0 = np.clip(struct0, 0, 1)
    return (1/TS*IS*struct) / struct0


def pt_yield_sugar(TS, IS, liqdens, conc, struct0, ctype=None):
    Mstruct, Msugar, Mdeg = _mol_weights(ctype)
    # enforce boundaries
    conc = np.clip(conc, 0, None)
    struct0 = np.clip(struct0, 0, 1)
    return (conc*(1-IS)/1000./Msugar/TS/liqdens) / (struct0/Mstruct)


def pt_yield_deg(TS, IS, TSr, liqdens, conc_l, conc_f, struct0, ctype,
                proportioned=False, syields=None, sfeeds=None):
    # mol weights are used to account for addition of water during hydrolysis
    # and subsequent degradation
    Mstruct, Msugar, Mdeg = _mol_weights(ctype)
    # enforce boundaries
    conc_l = np.clip(conc_l, 0, None)
    conc_f = np.clip(conc_f, 0, None)
    struct0 = np.clip(struct0, 0, 1)

    if proportioned:
        # enforce boundaries
        syields = np.clip(syields, 0, 1)
        sfeeds = np.clip(sfeeds, 0, 1)
        P = (1-syields[0])*sfeeds[0] / np.sum((1-syields)*sfeeds, axis=0)
    else:
        P = 1
    
    return P * (conc_l/TS*(1-IS)/liqdens + conc_f*(1/TSr-1/TS))/1000./Mdeg \
               / (struct0/Mstruct)


class Pretreater(Reactor):
    """
    A dilute-acid pretreatment group of variables and results. The object
    can calculate component yields of different structohydrates. Attributes
    may be set partially or fully, with calculations only executed if
    enough attributes are set for each calculated attribute. Components
    are set via dictionaries with keys for each structohydrate group,
    i.e. 'xylan'. Support for lignin is given by passing them under the
    guise of sugars (mono sugars is used for solubulized lignin).
    
    Assignable Attributes
    ---------------------
        Tr : float, (C)
            Reactor temperature. Estimate this from steam pressure.
        Tf : float, (C)
            Flash tank temperature.
        TS : float
            Flash tank slurry total solids frac.       
        IS : float
            Flash tank slurry insoluble solids frac (whole basis).
        TSr : float
            Total solids frac in the reactor; estimated by
            thermodynamic calculation method if not given.       
        cps : float, (kJ/kg-K)
            Specific heat of solids (assumed that of corn stover if not given).     
        liqdens : float, (g/mL)
            Density of slurry liquor.        
        struct0 : dict
            Feedstock structural structohydrate or lignin fractions in insol solids.       
        struct : dict
            Slurry structural structohydrate/lignin fractions in insol solids.   
        sug_total : dict, (g/L)
            Slurry liquor total sugar conccentrations.  
        sug_mono : dict, (g/L)
            Slurry liquor mono sugar or soluble lignin conccentrations.
        sug_oligo : dict, (g/L)
            Slurry liquor oligomeric sugar conccentrations (g/mL);
            calculated if not given.  
        deg_products : dict
            Assignment of degradation products to their sources, e.g.
            {'xylan':'furfural', 'arabinan':'furfural', 'glucan':'HMF'}
        deg : dict, (g/L)
            Degradation product concentrations in slurry liquor.
            Dictionary keys should match those specified in deg_products.            
        deg_f : dict, (g/L)
            Degradation product concentrations in flash tank condensate.
            Dictionary keys should match those specified in deg_products.            
        ctype : dict of strings
            Types of structohydrates: `pentose`, `hexose`, `lignin`;
            automatically set if struct0 dictionary keys are names of typical
            structohydrate families (xylan, arabinan, glucan, lignin, etc.).                
        proportioned : bool
            Option to set degradation product yields proportionately
            based on struct structohydrate yields, e.g. furfural is produced from
            degradation of both xylose and arabinose; therefore, furfural yield
            can be proportionately assigned to the source (xylan or arabinan)
            and thus mole closure is conserved.
                
    Calculated Attributes
    ---------------------
        y_struct : dict
            Yields of structural carbohydrates (molar basis) & lignin (mass basis). 
        y_struct_overall : float
            Total yield of all structural carbohydrates (excludes lignin).
        y_mono : dict
            Yields of monomeric sugars (molar basis) / solubulized lignin
            (mass basis).        
        y_mono_overall : float
            Total yield of all monomeric sugars (excludes lignin).
        y_oligo : dict
            Yields of oligomeric sugars (molar basis).        
        y_oligo_overall : float
            Total yield of all oligomeric sugars.
        y_deg : dict
            Yields of degradation products (molar basis).        
        y_deg_overall : float
            Total yield of all degradation products.
        mc : dict
            Mole closure for each component.
        omc : float
            Overall mole closure.
        
    Methods
    -------
        calcyields
            Run all yield calculations.
        plot
            Bar plot of the yields.
        export
            Export all attributes to excel.
    """    
    
    def __init__(self, **kwargs):
        # defaults
        self.proportioned = False
        self.cps = 1.463  # kJ/kg-K for corn stover
        # populate with kwarg settings
        allowed_keys = set([
                            'Tr', 
                            'Tf', 
                            'TS', 
                            'IS', 
                            'TSr', 
                            'liqdens',
                            'struct0',
                            'struct', 
                            'sug_total',
                            'sug_mono',
                            'sug_oligo', 
                            'deg', 
                            'deg_f', 
                            'ctype',
                            'proportioned', 
                            'cps', 
                            'deg_products',
                            'weighbelt',
                            'water1',
                            'water2',
                            'squeezate',
                            'steam',
                            'vent',
                            'slurry_flow',
                            'condensate',
                            ])
        for k, v in kwargs.items():
            if k in allowed_keys:
                setattr(self, k, v)
            else:
                raise TypeError('{} is an invalid keyword argument'.format(k))

    def calcstruct(self):
        self.y_struct, self.y_struct_overall = self.calcy(
                                                pt_yield_struct, 
                                                self.struct.keys(),
                                                TS=self.TS,
                                                IS=self.IS, 
                                                struct=self.struct, 
                                                struct0=self.struct0
                                                )

    def calcmono(self):
        self.y_mono, self.y_mono_overall = self.calcy(
                                                pt_yield_sugar, 
                                                self.sug_mono.keys(),
                                                TS=self.TS,
                                                IS=self.IS,
                                                liqdens=self.liqdens,
                                                conc=self.sug_mono,
                                                struct0=self.struct0,
                                                ctype=self.ctype,
                                                )


    def calcoligo(self):
        # calcualte oligos if they aren't there already
        if not hasattr(self, 'sug_oligo'):
            self.sug_oligo = {}
            for k in self.sug_total.keys():
                self.sug_oligo[k] = self.sug_total[k] - self.sug_mono[k]
        self.y_oligo, self.y_oligo_overall = self.calcy(
                                                pt_yield_sugar, 
                                                self.sug_oligo.keys(),
                                                TS=self.TS,
                                                IS=self.IS,
                                                liqdens=self.liqdens,
                                                conc=self.sug_oligo,
                                                struct0=self.struct0,
                                                ctype=self.ctype,
                                                )

    def calctotal(self):
        self.y_total, self.y_total_overall = self.calcy(
                                                pt_yield_sugar, 
                                                self.sug_total.keys(),
                                                TS=self.TS,
                                                IS=self.IS,
                                                liqdens=self.liqdens,
                                                conc=self.sug_total,
                                                struct0=self.struct0,
                                                ctype=self.ctype,
                                                )

    def calcdeg(self, TS, IS, liqdens, y_struct):
        # calculate degradation yields
        # calcualte TSr if not provided
        if not hasattr(self, 'TSr'):
            self.TSr = rx_solids(TS, self.Tr, self.Tf, self.cps)
        self.y_deg = {}
        syields = None
        sfeeds = None
        num, den = 0.0, 0.0
        for k in self.deg_products.keys():
            if self.proportioned:
                # build the proportioning array
                # component of interest is the first in the array
                syields = np.asarray([y_struct[k]])
                sfeeds = np.asarray([self.struct0[k]])
                # complete the array with everything else
                kset = set(y_struct.keys()) - set([k])
                for kk in kset:
                    if self.ctype[kk] == self.ctype[k]:  # only use same ctypes as the current deg_product
                        syields = np.vstack((syields, y_struct[kk]))
                        sfeeds = np.vstack((sfeeds, self.struct0[kk]))
            y_deg = pt_yield_deg(TS,
                                 IS, 
                                 self.TSr,
                                 liqdens, 
                                 self.deg[self.deg_products[k]],
                                 self.deg_f[self.deg_products[k]],
                                 self.struct0[k],
                                 self.ctype[k],
                                 self.proportioned,
                                 syields,
                                 sfeeds)
            # print(np.shape(syields))
            # print(np.shape(sfeeds))
            # note, `overall yield` does not consider lignin
            self.y_deg[k] = y_deg
            if self.ctype[k] != 'lignin':
                num, den = _molarizer(self, num, den, y_deg, k)
        self.y_deg_overall = num / den

    def TS_flow_correction(self):
        """
        Calculate what the total solids mass flow should be out of the reactor
        by accounting for water taken up during hydrolysis.
        (assumes total solids mass flow in is 1)
        """
        mols = 0
        for k in self.struct0.keys():
            if self.ctype[k] != 'lignin':
                M, _, _ = _mol_weights(self.ctype[k])
                mols += self.struct0[k] / M
        self.mTSs = 1 + (1-self.y_struct_overall) / self.TS * self.IS * 18 * mols

    def calcyields(self):
        if not hasattr(self, 'ctype'):
            self.setctype()
            
        if hasattr(self, 'struct0') and hasattr(self, 'struct') \
        and hasattr(self, 'TS') and hasattr(self, 'IS'):
            self.calcstruct()
            
        if hasattr(self, 'struct0') and hasattr(self, 'liqdens') \
        and hasattr(self, 'TS') and hasattr(self, 'IS') \
        and hasattr(self, 'sug_mono') and hasattr(self, 'ctype'):
            self.calcmono()
            
        if hasattr(self, 'struct0') and hasattr(self, 'liqdens') \
        and hasattr(self, 'TS') and hasattr(self, 'IS') \
        and ((hasattr(self, 'sug_mono') and hasattr(self, 'sug_total')) \
             or hasattr(self, 'sug_oligo')) \
        and hasattr(self, 'ctype'):
            self.calcoligo()
            
        if hasattr(self, 'struct0') and hasattr(self, 'liqdens') \
        and hasattr(self, 'TS') and hasattr(self, 'IS') \
        and hasattr(self, 'sug_total') and hasattr(self, 'ctype'):
            self.calctotal()
            
        if hasattr(self, 'struct0') and hasattr(self, 'liqdens') \
        and hasattr(self, 'TS') and hasattr(self, 'IS') \
        and hasattr(self, 'deg') and hasattr(self, 'deg_f') \
        and (hasattr(self, 'T1') or hasattr(self, 'TSr')) \
        and hasattr(self, 'ctype') and hasattr(self, 'deg_products'):
            self.calcdeg(self.TS, self.IS, self.liqdens, self.y_struct)
            
        if hasattr(self, 'y_struct') and hasattr(self, 'y_total') \
        and hasattr(self, 'y_deg'):
            self.calcmoleclose()
            
        self.TS_flow_correction()        
    
    def calcabsyield(self):
        self.y_mono_abs, self.y_mono_abs_overall = self.calcsug_abs(self.y_mono)
        self.y_oligo_abs, self.y_oligo_abs_overall = self.calcsug_abs(self.y_oligo)
        self.y_total_abs, self.y_total_abs_overall = self.calcsug_abs(self.y_total)
    
    def calcmassclose(self):
        massin = 0
        if hasattr(self, 'weighbelt'):
            massin += self.weighbelt
        if hasattr(self, 'water1'):
            massin += self.water1
        if hasattr(self, 'water2'):
            massin += self.water2
        if hasattr(self, 'chem'):
            massin += self.chem
        if hasattr(self, 'squeezate'):
            massin -= self.squeezate
        if hasattr(self, 'steam'):
            massin += self.steam
        massout = 0
        if hasattr(self, 'vent'):
            massout += self.vent
        if hasattr(self, 'slurry_flow'):
            massout += self.slurry_flow
        if hasattr(self, 'condensate'):
            massout += self.condensate
        
        if massin != 0:
            self.massclose = massout / massin
        else:
            print('Massclosure Error: no inlet flows provided.')


    


# =============================================================================
# Enzymatic Hydrolysis Yield Calculations
# =============================================================================

def eh_yield_struct(IS0, IS, struct0, struct, liqdens0=1, oligoconc0=0, ctype=None,
                    TS0=1, TS=1, flavor='proper'):
    """
    Leave TS0 and TS equal to 1 for true mass conservation calculations.
    TS0/TS can be used to apply the `total solids conservation` method when
    t0 data for EH, following dilution of the pretreated material, is not
    available (use post-pretreatment analytical values instead).
    """
    Mstruct, Msugar, Mdeg = _mol_weights(ctype)
    # enforce boundaries
    struct0 = np.clip(struct0, 0, 1)
    struct = np.clip(struct, 0, 1)
    oligoconc0 = np.clip(oligoconc0, 0, None)
    if flavor=='proper':
        return (IS/TS*struct/Mstruct) \
                / (IS0/TS0*struct0/Mstruct + (1-IS0)/TS0/liqdens0*oligoconc0/1000./Msugar)
    if flavor=='oldschool':  # ignore oligomers in the denominator
        return (IS/TS*struct/Mstruct) \
                / (IS0/TS0*struct0/Mstruct)

def eh_yield_sugar(IS0, IS, liqdens0, liqdens, conc0, conc,
                   struct0, oligoconc0=0, ctype=None, TS0=1, TS=1, flavor='proper'):
    """
    Leave TS0 and TS equal to 1 for true mass conservation calculations.
    TS0/TS can be used to apply the `total solids conservation` method when
    t0 data for EH, following dilution of the pretreated material, is not
    available (use post-pretreatment analytical values instead).
    """
    Mstruct, Msugar, Mdeg = _mol_weights(ctype)
    # enforce boundaries
    conc0 = np.clip(conc0, 0, None)
    conc = np.clip(conc, 0, None)
    oligoconc0 = np.clip(oligoconc0, 0, None)
    struct0 = np.clip(struct0, 0, 1)
    if flavor=='proper':
        return ((1-IS)/TS/liqdens*conc/1000./Msugar - (1-IS0)/TS0/liqdens0*conc0/1000./Msugar) \
                / (IS0/TS0*struct0/Mstruct + (1-IS0)/TS0/liqdens0*oligoconc0/1000./Msugar)
    if flavor=='oldschool':  # ignore oligomers in the denominator
        return ((1-IS)/TS/liqdens*conc/1000./Msugar - (1-IS0)/TS0/liqdens0*conc0/1000./Msugar) \
                / (IS0/TS0*struct0/Mstruct)

class Enzymer(Reactor):
    def __init__(self, **kwargs):
        # defaults
        self.method = 'masscons'
        self.flavor = 'proper'  # change to `oldschool` to ignore oligomers (inflate results)
        # populate with kwarg settings
        allowed_keys = set([
                            'TS0', 
                            'IS0',
                            'TS',
                            'IS',
                            'liqdens0',
                            'liqdens',
                            'struct0',
                            'struct', 
                            'sug_total0',
                            'sug_total',
                            'sug_mono0',
                            'sug_mono',
                            'sug_oligo0',
                            'sug_oligo',
                            'ctype',
                            'method',
                            'flavor',
                            ])
        for k, v in kwargs.items():
            if k in allowed_keys:
                setattr(self, k, v)
            else:
                raise TypeError('{} is an invalid keyword argument'.format(k))
        if self.method == 'masscons':
            # nullify TS effects for the formal mass conservation method
            self.TS0 = 1
            self.TS = 1

    def calcstruct(self):
        self.y_struct, self.y_struct_overall = self.calcy(
                                                eh_yield_struct, 
                                                self.struct1.keys(),
                                                IS0=self.IS0,
                                                IS=self.IS, 
                                                struct0=self.struct0, 
                                                struct=self.struct,
                                                liqdens0=self.liqdens0, 
                                                oligoconc0=self.sug_oligo0,
                                                ctype=self.ctype,
                                                TS0=self.TS0,
                                                TS=self.TS,
                                                flavor=self.flavor,
                                                )

    def calcmono(self):
        # calcualte oligos if they aren't there already
        if not hasattr(self, 'sug_oligo0'):
            self.sug_oligo0 = {}
            for k in self.sug_total0.keys():
                self.sug_oligo0[k] = self.sug_total0[k] - self.sug_mono0[k]
        self.y_mono, self.y_mono_overall = self.calcy(
                                                eh_yield_sugar, 
                                                self.sug_mono0.keys(),
                                                IS0=self.IS0,
                                                IS=self.IS, 
                                                liqdens0=self.liqdens0, 
                                                liqdens=self.liqdens, 
                                                conc0=self.sug_mono0, 
                                                conc=self.sug_mono,
                                                struct0=self.struct0, 
                                                oligoconc0=self.sug_oligo0,
                                                ctype=self.ctype,
                                                TS0=self.TS0,
                                                TS=self.TS,
                                                flavor=self.flavor,
                                                )

    def calcoligo(self):
        # calcualte oligos if they aren't there already
        if (not hasattr(self, 'sug_oligo1') or not hasattr(self, 'sug_oligo2')):
            self.sug_oligo0 = {}
            self.sug_oligo = {}
            for k in self.sug_total0.keys():
                self.sug_oligo0[k] = self.sug_total0[k] - self.sug_mono0[k]
                self.sug_oligo[k] = self.sug_total[k] - self.sug_mono[k]
        self.y_oligo, self.y_oligo_overall = self.calcy(
                                                eh_yield_sugar, 
                                                self.sug_oligo0.keys(),
                                                IS0=self.IS0,
                                                IS=self.IS, 
                                                liqdens0=self.liqdens0, 
                                                liqdens=self.liqdens, 
                                                conc0=self.sug_oligo0, 
                                                conc=self.sug_oligo,
                                                struct0=self.struct0, 
                                                oligoconc0=self.sug_oligo0,
                                                ctype=self.ctype,
                                                TS0=self.TS0,
                                                TS=self.TS,
                                                flavor=self.flavor,
                                                )

    def calctotal(self):
        self.y_total, self.y_total_overall = self.calcy(
                                                eh_yield_sugar, 
                                                self.sug_total0.keys(),
                                                IS0=self.IS0,
                                                IS=self.IS, 
                                                liqdens0=self.liqdens0, 
                                                liqdens=self.liqdens, 
                                                conc0=self.sug_total0, 
                                                conc=self.sug_total,
                                                struct0=self.struct0, 
                                                oligoconc0=self.sug_oligo0,
                                                ctype=self.ctype,
                                                TS0=self.TS0,
                                                TS=self.TS,
                                                flavor=self.flavor,
                                                )

    def calcyields(self):
        if not hasattr(self, 'ctype'):
            self.setctype()
        
        if hasattr(self, 'IS0') and hasattr(self, 'IS') \
        and hasattr(self, 'TS0') and hasattr(self, 'TS') \
        and hasattr(self, 'struct0'):
            if hasattr(self, 'struct'):
                self.calcstruct()

            if hasattr(self, 'liqdens0') and hasattr(self, 'liqdens'):
                if hasattr(self, 'sug_mono0') and hasattr(self, 'sug_mono'):
                    self.calcmono()
    
                if ((hasattr(self, 'sug_total0') and hasattr(self, 'sug_total')) \
                or (hasattr(self, 'sug_oligo0') and hasattr(self, 'sug_oligo'))):
                    self.calcoligo()
            
                if hasattr(self, 'sug_total0') and hasattr(self, 'sug_total'):
                    self.calctotal()
        if hasattr(self, 'y_struct'):
            self.calcmoleclose()



# =============================================================================
# Overall Yields Over Pretreatment and EH
# =============================================================================
class PTEHproc(Pretreater):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def calcstructPT(self):
        self.y_structPT, _ = self.calcy(
                                                pt_yield_struct, 
                                                self.structPT.keys(),
                                                TS=self.TSPT,
                                                IS=self.ISPT, 
                                                struct=self.structPT, 
                                                struct0=self.struct0
                                                )

    def calcyields(self):
        if not hasattr(self, 'ctype'):
            self.setctype()
            
        if hasattr(self, 'struct0') and hasattr(self, 'structPT') \
        and hasattr(self, 'TSPT') and hasattr(self, 'ISPT'):
            self.calcstructPT()
        if hasattr(self, 'struct0') and hasattr(self, 'struct') \
        and hasattr(self, 'TS') and hasattr(self, 'IS'):
            self.calcstruct()
            
        if hasattr(self, 'struct0') and hasattr(self, 'liqdens') \
        and hasattr(self, 'TS') and hasattr(self, 'IS') \
        and hasattr(self, 'sug_mono') and hasattr(self, 'ctype'):
            self.calcmono()
            
        if hasattr(self, 'struct0') and hasattr(self, 'liqdens') \
        and hasattr(self, 'TS') and hasattr(self, 'IS') \
        and ((hasattr(self, 'sug_mono') and hasattr(self, 'sug_total')) \
             or hasattr(self, 'sug_oligo')) \
        and hasattr(self, 'ctype'):
            self.calcoligo()
            
        if hasattr(self, 'struct0') and hasattr(self, 'liqdens') \
        and hasattr(self, 'TS') and hasattr(self, 'IS') \
        and hasattr(self, 'sug_total') and hasattr(self, 'ctype'):
            self.calctotal()
            
        if hasattr(self, 'struct0') and hasattr(self, 'liqdensPT') \
        and hasattr(self, 'TSPT') and hasattr(self, 'ISPT') \
        and hasattr(self, 'deg') and hasattr(self, 'deg_f') \
        and (hasattr(self, 'Tr') or hasattr(self, 'TSr')) \
        and hasattr(self, 'ctype') and hasattr(self, 'deg_products'):
            self.calcdeg(self.TSPT, self.ISPT, self.liqdensPT, self.y_structPT)
            
        if hasattr(self, 'y_struct') and hasattr(self, 'y_total') \
        and hasattr(self, 'y_deg'):
            self.calcmoleclose()
            

# =============================================================================
# Pretreatment Reactor Operating Conditions Planning
# =============================================================================




# =============================================================================
# Testing
# =============================================================================
if __name__ == '__main__':
    y = Pretreater(Tr=160, Tf=94, TS=0.35, IS=0.20, liqdens=1.05,\
                 struct0={'x':0.29, 'g':0.33, 'ara':0.10, 'lig':0.3},
                 struct={'x':0.05, 'g':0.28, 'ara':0.03, 'lig':0.45},
                 sug_total={'x':80, 'g':20, 'ara':15},
                 sug_mono={'x':60, 'g':18, 'ara':10, 'lig':15},
                 deg_products={'x':'f', 'g':'hmf', 'ara':'f'},
                 deg={'f':5, 'hmf':1},
                 deg_f={'f':15, 'hmf':2},
                 # ctype={'x':'pentose', 'g':'hexose', 'ara':'pentose'},
                 proportioned=True,
                 )
                 
    y.calcyields()
    # y.setctype()
    # y.calcstruct()
    # y.export('test.xlsx')
            
            


