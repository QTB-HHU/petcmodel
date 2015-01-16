# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 11:42:21 2014

Calculate the steady state of the system when state transition are switched off


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
from mpl_toolkits.mplot3d import Axes3D

import petcModel
import lightProtocol
import parametersPETC
import simulate
from misc import pH, pHinv

p = parametersPETC.ParametersPETC()
# ---------------------------------------- #
# switch off the state transitions #
p.staticAntI = 0
p.staticAntII = 0
p.kStt7 = 0
p.kPph1 = 0

m = petcModel.PETCModel(p)
s = simulate.Sim(m)

PFDrange = np.linspace(75,1500,20)
STrange = np.linspace(0,1,20)

Ys = np.zeros([len(STrange), len(PFDrange), 8])

# dark adapted state. Not important
y0=np.array([p.PQtot, 0.0202, 5.000, 0.0000, 0.0000, 0.0001, 0.9, 0.0000])

for i in range(len(STrange)):
	y0[6] = STrange[i]
	print('Teraz liczymy dla ' + str(i)) 
	Y = s.steadyStateLightScan(PFDrange,y0)
	Ys[i,:] = Y

import pickle
ss = {"STrange": STrange, "PFDrange": PFDrange, "Ys": Ys}
output = open('steadyStateAnalysisFixedST.pkl', 'wb')
pickle.dump(ss, output, 2)
output.close()



#================================================================= #
# 3 D plotting routine to obtain figure as in Ebenhoeh et al. 2014 #
'''
import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

file = open('steadyStateAnalysisFixedST_corrected.pkl', 'rb')
data = pickle.load(file)


ST = data['STrange']
PFD = data['PFDrange']
X,Y = np.meshgrid(PFD, ST)

Yss = np.zeros([len(ST),20])

for i in range(len(ST)):
        for j in range(len(PFD)):
         Yss[i,j] = 1 - data['Ys'][i][j][0] / 17.5

y_formatter = matplotlib.ticker.ScalarFormatter(useOffset = False)
ax.yaxis.set_major_formatter(y_formatter)

cm = matplotlib.cm.get_cmap('jet')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X,Y,Yss, rstride=1, cstride=1, cmap=cm,
                       linewidth=1, antialiased=True)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

ax.set_zlim3d(0, 1)
fig.colorbar(surf) #, shrink=0.5, aspect=5)
plt.title('steady state of reduced plastoquinon pool')
plt.xlabel('PFD')
plt.ylabel('PSII cross section')
plt.show()
'''