#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 16:53:42 2014

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
import matplotlib.pyplot as plt

import petcModel
import lightProtocol
import parametersPETC
import simulate
from misc import pH, pHinv
import petcResults

# ================================================================ #
p = parametersPETC.ParametersPETC()
m = petcModel.PETCModel(p)
s = simulate.Sim(m)

# declare desired light protocol
l = lightProtocol.LightProtocol({'protocol':'PAM_dld','PFD':0.001, 'Tflash':0.7, 'dtFlash':120, 'Ton':300, 'Toff':1600})

# load dark adaptation or start with this initial condition
y0=np.array([p.PQtot, 0.0202, 5.000, 0.0000, 0.0000, 0.0001, 0.9, 0.0000, 0.0000])

# declare lenght of the simulation [in s]
Tmax = 3500 

t0 = 0
t1 = l.Tfl
cnt = 0

# =========================================== #
# Anoxia experiment simulation in three parts #

# ======= 1 ======== #
# first aerobic part #
while t1 < l.Ton:
    time = np.linspace(t0,t1,101)
    Y = s.timeCourse(l,time,y0)
    y0 = Y[-1,:]
    cnt += 1
    t0 = t1
    if cnt % 2 == 0:
        t1 = t0 + l.Tfl
    else:
        t1 = t0 + (l.dt - l.Tfl)
        
# overal results of this model put pass to the BigR array	     
BigR = s.results
cnt =0

# ======= 2 ======== #
# next anerobic part #

p.O2ext = 0
m = petcModel.PETCModel(p)
s2 = simulate.Sim(m)
     
while t1 >= l.Ton and t1 < l.Toff:
    time = np.linspace(t0,t1,101)
    
    Y = s2.timeCourse(l,time,y0)
    y0 = Y[-1,:]
    cnt += 1
    t0 = t1
    if cnt % 2 == 0: #parzyste
        t1 = t0 + l.Tfl
    else:
        t1 = t0 + (l.dt - l.Tfl)

BigR = np.hstack([BigR, s2.results])
cnt = 0

# ========= 3 ========= #
# fback to aerobic part #

p.O2ext = 8
p.kNDH = 0
m = petcModel.PETCModel(p)
s3 = simulate.Sim(m)

while t1 >= l.Toff and t1 < Tmax:
    time = np.linspace(t0,t1,101)
    Y = s3.timeCourse(l,time,y0)
    y0 = Y[-1,:]
    cnt += 1
    t0 = t1
    if cnt % 2 == 0: #parzyste
        t1 = t0 + l.Tfl
    else:
        t1 = t0 + (l.dt - l.Tfl)

BigR = np.hstack([BigR, s3.results])

s.results = BigR 

# ============================ #
# Simulate fluorescence traces #

res = petcResults.PETCResults(s)

res.plotFluo()
plt.show()