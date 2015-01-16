# File parametersPETC.py
# -*- coding: utf-8 -*-
"""

Copyright (C) 2014-2015  Anna Matuszyńska, Oliver Ebenhöh

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (license.txt).  If not, see <http://www.gnu.org/licenses/>.
"""

from numpy import log, exp

class ParametersPETC:

    defaultparameterset = {
        # pool sizes
        'PSIItot': 2.5, # [mmol/molChl] total concentration of PSII
        'PSItot': 2.5,
	'PQtot': 17.5, # [mmol/molChl]
	'PCtot': 4., # Bohme1987, but other sources give different values - seems to depend greatly on organism and conditions
        'Fdtot': 5., # Bohme1987
        'Ctot': 2.5, #source unclear (Schoettler says 0.4...?, but plausible to assume that complexes (PSII,PSI,b6f) have approx. same abundance)
        'NADPtot': 25., # estimate from ~ 0.8 mM, Heineke1991
	'APtot': 60., # [mmol/molChl] Bionumbers ~2.55mM (=81mmol/molChl) (FIXME: Soma had 50)

	
        # parameters associated with photosystem II
        'kH': 0.,
        'kH0': 5.e8, # base quenching' after calculation with Giovanni
        'kF': 6.25e7, # fluorescence 16ns
        'k1': 5.e9, # excitation of Pheo / charge separation 200ps
        'k1rev': 1.e10,
        'k2': 5.e9, # original 5e9 (charge separation limiting step ~ 200ps) - made this faster for higher Fs fluorescence

        # parameters associated with photosystem I
        'kStt7': 0.0035, # [s-1] fitted to the FM dynamics
        'kPph1': 0.0013, # [s-1] fitted to the FM dynamics
        'KM_ST': 0.2, # Switch point (half-activity of Stt7) for 20% PQ oxidised (80% reduced)
        'n_ST': 2., # Hill coefficient of 4 -> 1/(2.5^4)~1/40 activity at PQox=PQred
        'staticAntI': 0.2,
        'staticAntII': 0.0,

        # ATP and NADPH parameters
        'kATPsynth': 20., # taken from MATLAB
        'kATPcons': 10., # taken from MATLAB
        'kATPimport': 0., # TODO possibility for ATP import at night - NOT YET IMPLEMENTED!
        'ATPcyt': 0.5, # only relative levels are relevant (normalised to 1) to set equilibrium
        'Pi_mol': 0.01,
        'DeltaG0_ATP': 30.6,  # 30.6kJ/mol / RT
        'HPR': 14./3.,
        'kNADPHimport': 0., # TODO possibility for NADPH import - NOT YET IMPLEMENTED!
        'kNADPHcons': 15., # taken from MATLAB
        'NADPHcyt': 0.5, # only relatice levels

        # global conversion factor of PFD to excitation rate
        'cPFD': 4., # [m^2/mmol PSII] 

        # pH and protons
        'pHstroma': 7.8,
        'kLeak': 0.010, # [1/s] leakage rate -- inconsistency with Kathrine
        'bH': 100., # proton buffer: ratio total / free protons

        # rate constants
        'kPQred': 250., # [1/(s*(mmol/molChl))]
        'kCytb6f': 2.5, # a rough estimate: transfer PQ->cytf should be ~10ms
        'kPTOX': .01, # ~ 5 electrons / seconds. This gives a bit more (~20)
        'kPCox': 2500., # a rough estimate: half life of PC->P700 should be ~0.2ms
        'kFdred': 2.5e5, # a rough estimate: half life of PC->P700 should be ~2micro-s
        'kcatFNR': 500., # Carrillo2003 (kcat~500 1/s)
        'kcyc': 1.,

        'O2ext': 8., # corresponds to 250 microM, corr. to 20%
        'kNDH': .002, # re-introduce e- into PQ pool. Only positive for anaerobic (reducing) condition
        'kNh': 0.05,
        'kNr': 0.004,
        'NPQsw': 5.8,
        'nH': 5.,

        'EFNR': 3., # Bohme1987
        'KM_FNR_F': 1.56, # corresponds to 0.05 mM (Aliverti1990)
        'KM_FNR_N': 0.22, # corresponds to 0.007 mM (Shin1971, Aliverti2004)

        # standard redox potentials (at pH=0) in V
        'E0_QA': -0.140,
        'E0_PQ': 0.354,
        'E0_cytf': 0.350,
        'E0_PC': 0.380,
        'E0_P700': 0.480,
        'E0_FA': -0.550,
        'E0_Fd': -0.430,
        'E0_NADP': -0.113,

        # physical constants
        'F': 96.485, # Faraday constant
        'R': 8.3e-3, # universal gas constant
        'T': 298., # Temperature in K - for now assumed to be constant at 25 C
    }


    def __init__(self, pars = {}):
        mypars = pars.copy()
        for k in ParametersPETC.defaultparameterset.keys():
            mypars.setdefault(k,ParametersPETC.defaultparameterset[k])

        for k in mypars.keys():
            setattr(self,k,mypars[k])

        self.setCompositeParameters()

        setattr(self,'KeqPQred',self.Keq_PQred())
        setattr(self,'KeqCyc', self.Keq_cyc())
        setattr(self,'KeqCytfPC', self.Keq_cytfPC())
        setattr(self,'KeqFAFd', self.Keq_FAFd())
        setattr(self,'KeqPCP700', self.Keq_PCP700())
        setattr(self,'KeqNDH', self.Keq_NDH())
        setattr(self,'KeqFNR', self.Keq_FNR())


    def setCompositeParameters(self):
        setattr(self, 'RT', self.R * self.T)
        setattr(self, 'dG_pH', log(10)*self.RT)
        setattr(self, 'Hstroma', 3.2e4*10**(-self.pHstroma)) # proton concentration in stroma
        setattr(self, 'kProtonation', 4e-3 / self.Hstroma) # [1/s] converted from 4 * 10^-6 [1/ms] protonation of LHCs (L), depends on pH value in lumen


    def Keq_PQred(self):
        DG1 = -self.E0_QA * self.F
        DG2 = -2 * self.E0_PQ * self.F
        DG = -2 * DG1 + DG2 + 2 * self.pHstroma * self.dG_pH
        K = exp(-DG/self.RT)
        return K

    def Keq_cyc(self):
        DG1 = -self.E0_Fd * self.F
        DG2 = -2 * self.E0_PQ * self.F
        DG = -2 * DG1 + DG2 + 2 * self.dG_pH * self.pHstroma
        K = exp(-DG/self.RT)
        return K

    def Keq_cytfPC(self):
        DG1 = -self.E0_cytf * self.F
        DG2 = -self.E0_PC * self.F
        DG = -DG1 + DG2
        K = exp(-DG/self.RT)
        return K

    def Keq_FAFd(self):
        DG1 = -self.E0_FA * self.F
        DG2 = -self.E0_Fd * self.F
        DG = -DG1 + DG2
        K = exp(-DG/self.RT)
        return K

    def Keq_PCP700(self):
        DG1 = -self.E0_PC * self.F
        DG2 = -self.E0_P700 * self.F
        DG = -DG1 + DG2
        K = exp(-DG/self.RT)
        return K

    def Keq_NDH(self):
        DG1 = -2 * self.E0_NADP * self.F
        DG2 = -2 * self.E0_PQ * self.F
        DG = -DG1 + DG2 + self.dG_pH * self.pHstroma
        K = exp(-DG/self.RT)
        return K

    def Keq_FNR(self):
        DG1 = -self.E0_Fd * self.F
        DG2 = -2 * self.E0_NADP * self.F
        DG = -2 * DG1 + DG2 + self.dG_pH * self.pHstroma
        K = exp(-DG/self.RT)
        return K
