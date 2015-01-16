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

p = parametersPETC.ParametersPETC()
m = petcModel.PETCModel(p)
s = simulate.Sim(m)

l = lightProtocol.LightProtocol({'protocol':'PAM_dld','PFD':100, 'Tflash':0.7, 'dtFlash':90, 'Ton':270, 'Toff':900})

y0=np.array([p.PQtot, 0.0202, 5.000, 0.0000, 0.0000, 0.0001, 1, 0.0000])
Tmax = 2100

t0 = 0
t1 = l.Tfl
cnt = 0

while t0 < Tmax:
     time = np.linspace(t0,t1,101)
     Y = s.timeCourse(l,time,y0)
     y0 = Y[-1,:]
     cnt += 1
     t0 = t1
     if cnt % 2 == 0:
        t1 = t0 + l.Tfl
     else:
        t1 = t0 + (l.dt - l.Tfl)

import petcResults
res = petcResults.PETCResults(s)

# plot fluorescence trace #
res.plotFluo()
plt.show()