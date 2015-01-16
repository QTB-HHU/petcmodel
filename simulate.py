# File simulate.py
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

import numpy as np
from scipy.integrate import ode
from math import floor

import petcModel
import lightProtocol
import parametersPETC
from misc import pH, pHinv

class Sim:
    """
    class includes integration methods
    """
    def __init__(self, model):

        self.model = model

        def dydt(t, y, m, l):
            return m.petc_model(y,l.lightintensity(t))

        self.dydt = dydt
        
        self._successful = True
        
        self._monitor = True
        
        self._warnings = False
        
        self.clearResults()
        
    
    def successful(self):
        return self._successful

    def doesMonitor(self, setMonitor=None):
        if setMonitor != None:
            self._monitor = setMonitor
        return self._monitor
    
    def clearResults(self):
        self.results = []
    

    # integration, returns variables at time t
    def integrate(self, lightFn, t, y0, integrator='lsoda', minstep=1e-8, maxstep=0.1, nsteps=500, t0=0):

        step = maxstep

        numsteps = max(nsteps, 10*floor((t-t0)/step))
 
        while step >= minstep:
            r = ode(self.dydt).set_integrator(integrator,max_step=step,nsteps=numsteps)
            r.set_initial_value(y0,t0)
            r.set_f_params(self.model, lightFn)

            # suppress FORTRAN warnings
            if not self._warnings:
                r._integrator.iwork[2] = -1

            try:
                r.integrate(t)
                
                if r.successful():
                    break
                
            except petcModel.ModelError:
                print('caught error at ',step,'. Reducing step size')
                #step = step/10
            
            step = step/10
            numsteps = numsteps*10
            
            if self._warnings:           
                print('numsteps=',numsteps,', step=',step)
                print(r.t,r.y)
                print(self.model.rates(r.y,lightFn.lightintensity(r.t)))
                
        self._successful = r.successful()
        return r.y

    def timeCourse(self, lightFn, T, y0, integrator='lsoda', minstep=1e-8, maxstep=0.1, nsteps=500):
        ''' integration over time, different integrators possible, lsoda default
            returns: array of state variables
        '''
        self._successful = True
        
        Y = [y0]

        cnt = 1
        while cnt<len(T) and self.successful():
            Y.append(self.integrate(lightFn,T[cnt],Y[cnt-1],t0=T[cnt-1], 
                                    minstep=minstep, 
                                    maxstep=maxstep, 
                                    nsteps=nsteps,
                                    integrator=integrator))
            print(T[cnt])
            cnt+=1

        if self.doesMonitor() and self.successful():
            self.results.append({'t':T,'y':np.vstack(Y),'lfn':lightFn})
        
        return np.vstack(Y)


    def steadyStateLight(self, PFD, AbsTol = 1e-3, 
                         Tstep = 1,
                         maxstep = 1000,
                         y0=np.array([8.75, 0.0202, 5.000, 0.0000, 0.0000, 0.0001, 0.9, 0.0000, 0.0000])
                         ):
        
        l = lightProtocol.LightProtocol({'protocol':'const','PFD':PFD})
        
        T = 0
        cnt = 0
        Y0 = y0
        err = np.linalg.norm(y0,ord=2)
        
        while self.successful() and cnt < maxstep and err > AbsTol:
            Y = self.integrate(l,T+Tstep,Y0,t0=T)
            T += Tstep
            cnt += 1
            err = np.linalg.norm(Y-Y0,ord=2)
            print('T=',T,' err=',err)
            Y0 = Y
        
        return Y
        
            
                             
        
        
        
    def steadyStateLightScan(self, PFDrange, y0, t = 1000):
        ''' steady state at time t for given PFDs '''

        Ys = np.zeros((len(PFDrange),len(y0)))
        T = np.linspace(0,t,t*10) # time course method needs vector of time

        for cnt in range(len(PFDrange)):
            PFD = PFDrange[cnt]
            l = lightProtocol.LightProtocol({'protocol':'const','PFD':PFD})

            Ys[cnt,:] = self.timeCourse(l, T, y0)[-1,:]

        return Ys