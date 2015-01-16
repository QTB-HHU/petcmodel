#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Full model of the photosynthetic electron transport chain

petcModel defines methods to calculate reaction rates and set of eight ODEs
based on the model published by Ebenhoeh et al. in 2014


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


import numpy as np
from numpy import log,exp
from misc import pH, pHinv

import lightProtocol

__author__ = "Anna Matuszyńska"
__copyright__ = "Copyright 2014, Heinrich-Heine University Dusseldorf"
__credits__ = ["Anna Matuszynska", "Oliver Ebenhoeh"]
__maintainer__ = "Anna Matuszynska"
__email__ = "Anna.Matuszynska@uni-duesseldorf.de"
__status__ = "Development"

class PETCModel(object):

    def __init__(self, par):
        """ returns set of parameters """
        self.par = par
        
        print('You run simulations based on the model published in 2014 by Ebenhoeh et al.')

    # -------------------------------------------- #
    # methods required to calculate reaction rates #

    def ps2states(self, Pox, Q, L):
        """ QSSA, calculates the states of the photosystem II
            accepts values:
            Pox: oxidised fraction of the PQ pool (PQH2)
            Q: quencher
            L: light, int or array of the n x 1 dimension, that gives the light intensity

            returns:
            B: array of arrays with the states of PSII; rows: time, columns states: 1 and 3 excited
        """
        Pred = self.par.PQtot - Pox

        k2 = self.par.k2
        kF = self.par.kF
        kH = self.par.kH0 + self.par.kH * Q

        k3p = self.par.kPQred * Pox
        k3m = self.par.kPQred * Pred / self.par.KeqPQred;

        M = np.array([[-L-k3m, kH+kF,       k3p, 0],
                      [L,      -(kH+kF+k2), 0,   0],
                      [0,      0,           L,   -(kH+kF)],
                      [1,      1,           1,   1]])

        A = np.array([0, 0, 0, self.par.PSIItot])

        # B = [B0,B1*,B2,B3*]
        B = np.linalg.solve(M, A)

        return B

    def ps1states(self, C, F, L):
        """ QSSA calculates open state of PSI
        depends on reduction states of plastocyanin and ferredoxin
        C = [PC], F = [Fd] (ox. forms)
        accepts: light, y as an array of arrays
        returns: array of PSI open
        """

        #C = y[:,1]
        #F = y[:,2]

        A1 = self.par.PSItot / \
        (1 + L/(self.par.kFdred * F) + \
        (1 + (self.par.Fdtot - F)/(self.par.KeqFAFd * F)) * (C/(self.par.KeqPCP700 * (self.par.PCtot - C)) + \
        L/(self.par.kPCox * (self.par.PCtot-C)))
        )

        return A1


    # ---------------------------------- #
    # method to calculate cross sections #

    def crossSectionSimple(self, Ant):
        """ calculates the cross section of PSII """
        cs = self.par.staticAntII + (1 - self.par.staticAntII - self.par.staticAntI) * Ant
        return cs


    # --------------------------- #
    # 14 reaction rates seperatly #

    def vPS2(self,Pox,H,Q,LII):
        """ reaction rate constant for photochemistry """
        Q = self.quencher(Q,H)

        B = self.ps2states(Pox,Q,LII)

        v = self.par.k2 * B[1] / 2
        return v


    def vPS1(self, C,F,LI):
        """ reaction rate constant for open PSI """
        A = self.ps1states(C,F,LI)
        v = LI * A
        return v


    def vPTOX(self, Pox):
        """ calculates reaction rate of PTOX """
        v = (self.par.PQtot - Pox) * self.par.kPTOX * self.par.O2ext
        return v


    def vB6f(self, Pox,C,H):
        """ calculates reaction rate of cytb6f """
        ph = pH(H)
        Keq = self.Keq_cytb6f(ph)

        v = max(self.par.kCytb6f * ((self.par.PQtot - Pox) * C**2 - (Pox * (self.par.PCtot - C)**2)/Keq), -self.par.kCytb6f)
        return v


    def vNDH(self, Pox):
        """ calculates reaction rate of PQ reduction under absence of oxygen
            can be mediated by NADH reductase NDH """
        v = self.par.kNDH * Pox
        return v


    def vCyc(self, Pox, F):
        """ calculates reaction rate of cyclic electron flow
            considered as practically irreversible """
        v = self.par.kcyc * (((self.par.Fdtot - F)**2) * Pox)
        return v


    def vFNR(self, F,N):
        """ reaction rate mediated by FNR
            uses convenience kinetics """
        fdred = (self.par.Fdtot - F)/self.par.KM_FNR_F
        fdox = F/self.par.KM_FNR_F
        nadph = N/self.par.KM_FNR_N
        nadp = (self.par.NADPtot - N)/self.par.KM_FNR_N

        v = (self.par.EFNR * self.par.kcatFNR *
            ((fdred**2) * nadp - ((fdox**2) * nadph) / self.par.KeqFNR) /
            ((1+fdred+fdred**2) * (1+nadp) + (1+fdox+fdox**2) * (1+nadph) - 1))
        return v


    def vLeak(self, Hlf):
        """ rate of leak of protons through the membrane """
        v = self.par.kLeak * (Hlf - pHinv(self.par.pHstroma))
        return v


    def vSt12(self, Pox,Ant):
        """ reaction rate of state transitions from PSII to PSI
            Ant depending on module used corresponds to non-phosphorylated antennae
            or antennae associated with PSII
        """
        kKin = self.par.kStt7 * ( 1 / (1 + ((Pox /self.par.PQtot)/self.par.KM_ST)**self.par.n_ST))
        v = kKin * Ant
        return v


    def vSt21(self, Ant):
        """ reaction rate of state transitions from PSI to PSII """
        v = self.par.kPph1 * (1 - Ant)
        return v


    # ---------------------- #
    # rates of main products #

    def vATPsynthase(self, A, Hlf):
        ph = pH(Hlf)
        ADP = self.par.APtot - A
        v = self.par.kATPsynth * (ADP - A / self.Keq_ATP(ph))
        return v


    def vATPconsumption(self, A):
        v = self.par.kATPcons * A
        return v


    def vNADPHconsumption(self, N):
        v = self.par.kNADPHcons * N
        return v

    # ---------------------------------------------- #
    # model (parameter) dependent helper functions   #

    def Keq_ATP(self, pH):
        DG = self.par.DeltaG0_ATP - self.par.dG_pH * self.par.HPR * (self.par.pHstroma - pH)
        Keq = self.par.Pi_mol * exp(-DG/self.par.RT)
        return Keq


    def Keq_cytb6f(self, pH):
        DG1 = -2 * self.par.F * self.par.E0_PQ
        DG2 = -self.par.F * self.par.E0_PC

        DG = - (DG1 + 2*self.par.dG_pH * pH) + 2 * DG2 + 2*self.par.dG_pH * (self.par.pHstroma - pH)
        Keq = exp(-DG/self.par.RT)
        return Keq

    def quencher(self, Q, H):
        vNPQdyn = H**self.par.nH / (H ** self.par.nH + pHinv(self.par.NPQsw)**self.par.nH)
        v = (1-Q) * self.par.kNh * vNPQdyn - self.par.kNr * Q
        return v

    # ----------------------------------------------#
    # reaction rates all together

    def rates(self, y, PFD):
        """ Method calculating reaction rates.
            called as m.rates(y,PFD)
            Input: y - vector of state variables and PFD - photon flux density
        """

        P = y[0] # oxidised plastoquinone
        C = y[1] # plastocyan
        F = y[2] # ferrodoxin
        A = y[3] # concentration of ATP
        N = y[4] # concentration of NADPH
        H = y[5] # lumenal protons
        T = y[6] # non-phosphorylated antenna // associated with PSII
        Q = y[7] # quencher
        
        cs2 = self.crossSectionSimple(T)
        LII = cs2 * PFD
        LI = (1-cs2) * PFD # assumption: total cross section is always =1


        v = {
            'ps2': self.vPS2(P,H,Q,LII), 

            'ps1': self.vPS1(C,F,LI),

            'ptox': self.vPTOX(P),

            'cytb6f': self.vB6f(P,C,H),

            'ndh': self.vNDH(P),

            'fnr': self.vFNR(F,N),

            'cyc': self.vCyc(P,F),

            'ATPsynthase': self.vATPsynthase(A,H),

            'ATPconsumption': self.vATPconsumption(A),

            'NADPHconsumption': self.vNADPHconsumption(N),

            'leak': self.vLeak(H),

            'st12': self.vSt12(P,T),

            'st21': self.vSt21(T),

            'quencher': self.quencher(Q,H)

            }

        return v


    def petc_model(self, y, PFD):

            '''' Defining the system of NINE equations governing the evolution of the photosynthetic ETC
                Returns an array containing the value of y for each desired time in t, with the initial value `y0` in the first row.
                Updated on 28 Oct 2014
            '''
            # error kept as in npq_model as similiar problems expected
            if y[5] <= 0:
                raise ModelError("H went below zero")

            # Reaction rates, v is a dictionary
            v = self.rates(y, PFD)

            # Output from ODE function must be a COLUMN vector, with n rows
            n = len(y)      # 1: implies its a single ODE
            dydt = np.zeros((n,1))

            # dP/dt = -vPSII + vCytb6f - vcyc + vPTOX - vNDH
            dydt[0] = - v['ps2'] + v['cytb6f'] - v['cyc'] + v['ptox'] - v['ndh']
            # dC/dt = -2vCytb6f + vPSI
            dydt[1] = -2 * v['cytb6f'] + v['ps1']
            # dF/dt = -vPSI + 2vFNR + 2vcyc
            dydt[2] = - v['ps1'] + 2 * v['fnr'] + 2 * v['cyc']
            # dATP/dt = vATPsynthase - vATPconsumption + TODO: import
            dydt[3] = v['ATPsynthase'] - v['ATPconsumption']
            # dN/dt = vFNR - vNADPHconsumption + TODO: vNADPHimport
            dydt[4] = v['fnr'] - v['NADPHconsumption']
            # dH/dt = water_split + proton_pump - synthase_ATP - leak
            dydt[5] = (2 * v['ps2'] + 4 * v['cytb6f'] - self.par.HPR * v['ATPsynthase'] - v['leak']) / self.par.bH
            # dT/dt =  non-phosphorylated antennae // antenna associated with PSII
            dydt[6] = - v['st12'] + v['st21']
            # original MATLAB model with simpel quencher
            dydt[7] = v['quencher']
            
            return dydt

class ModelError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
