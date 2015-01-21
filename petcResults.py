# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 12:44:28 2014

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

import matplotlib.pyplot as plt
import numpy as np
from misc import pH
from simulate import Sim

class PETCResults(Sim):


    def __init__(self,s):

        self.model = s.model
        self.results = s.results
        self.par = s.model.par

    def plotph(self, r=None):
        if r == None:
            r = range(len(self.results))

        for i in r:
            plt.plot(self.results[i]['t'],pH(self.results[i]['y'][:,5]))


    def plotRel(self, r=None):
        """ plot results of integration
            six plots for PQ, PC, Fd, ATP, NADPH antennae
        """

        if r == None:
            r = range(len(self.results))

        for i in r:
            t = self.results[i]['t']
            y = self.results[i]['y']
            if i == 0:  # unique labels
                plt.plot(t, 1 - y[:,0]/self.par.PQtot, color='red', label = 'PQred') # PQred
                plt.plot(t, 1 - y[:,1]/self.par.PCtot, color='blue', label = 'PCred') # PC
                plt.plot(t, 1 - y[:,2]/self.par.Fdtot, color='yellow', label = 'Fdred') # Fd
                plt.plot(t, y[:,3]/self.par.APtot, color='magenta', label = 'ATP') # ATP
                plt.plot(t, y[:,4]/self.par.NADPtot, color='cyan', label = 'NADPH') # NADPH
                plt.plot(t, y[:,6], color='green', label = 'CSII') # antennae on photosystem II
            else:
                plt.plot(t, 1 - y[:,0]/self.par.PQtot, color='red')     # PQred
                plt.plot(t, 1 - y[:,1]/self.par.PCtot, color='blue')    # PC
                plt.plot(t, 1 - y[:,2]/self.par.Fdtot, color='yellow')  # Fd
                plt.plot(t, y[:,3]/self.par.APtot, color='magenta')     # ATP
                plt.plot(t, y[:,4]/self.par.NADPtot, color='cyan')      # NADPH
                plt.plot(t, y[:,6], color='green')                      # antennae on photosystem II

            plt.xlabel('time')
            plt.title('Temporal evolution of state variables')
            plt.legend(loc= 'best') 


    def plotV(self, v, r=None):

        if r == None:
            r = range(len(self.results))


        for i in r:
            t = self.results[i]['t']
            y = self.results[i]['y']
            l = self.results[i]['lfn']

            V = [self.model.rates(y[i],l.lightintensity(t[i]))[v] for i in range(len(t))]
            plt.plot(t,V)


    def fluo(self, r=None):
        '''
        :return: time vector, fluorescence, maximal fluorescence, steady state fluorescence
        '''

        if r == None:
            r = range(len(self.results))

        FM = []
        F0 = []
        F = []
        T = []

        for i in r:

            t = self.results[i]['t']
            y = self.results[i]['y']
            l = self.results[i]['lfn']
            P = y[:,0]
            anT = y[:,6] # NEVER USE T for antennae. it is reserved for time
            L = y[:,7]
            H = y[:,5]

            Q = self.model.quencher(L,H)
            cs2 = self.model.crossSectionSimple(anT)

            B = [self.model.ps2states(P[i], Q[i], cs2[i] * l.lightintensity(t[i])) for i in range(len(t))]

            Fzero = [cs2[i] * self.par.kF / (self.par.kF + self.par.kH0 + self.par.kH * Q[i] + self.par.k2)*x for x in [B[i][0] for i in range(len(t))]]
            Fm = [cs2[i] * self.par.kF / (self.par.kF + self.par.kH0 + self.par.kH * Q[i])*x for x in [B[i][2] for i in range(len(t))]]
            F_ = [x + y for x, y in zip(Fm, Fzero)]

            T = np.hstack([T,t])
            FM = np.hstack([FM, Fm])
            F0 = np.hstack([F0, Fzero])

        F = [x + y for x, y in zip(FM, F0)] # could hstack but was comparing two methods

        return T, F, FM, F0



    def plotFluo(self):
        '''
            plot fluorescence trace, uses function fluo to calculate FM and F0
        '''
        T, F, _, _ = self.fluo()
        plt.plot(T,F/max(F), 'r')
        plt.xlabel('time')
        plt.title('Fluorescence trace')
