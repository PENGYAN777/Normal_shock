#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 16:05:08 2022

@author: yan

plot (P,T) diagram. show contour of compressibility factor Z and 
fundamental derivative of gasdynamics Gamma
"""

import matplotlib
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import scipy.interpolate


# give fluid name
fluidname = "MM"

# update fluid
fluid = CP.AbstractState("HEOS", fluidname)

pc = fluid.keyed_output(CP.iP_critical)
Tc = fluid.keyed_output(CP.iT_critical)
Tmin =  CP.CoolProp.PropsSI("Ttriple",fluidname)
Tmax =  CP.CoolProp.PropsSI("Tmax",fluidname)
pmax = fluid.keyed_output(CP.iP_max)
pmin = fluid.keyed_output(CP.iP_min)
fillcolor = 'g'

# fig = plt.figure(figsize = (6,6))
# ax = fig.add_subplot(111)
fig = plt.figure(figsize=(6,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height]) 
lw = 3


# ----------------
# Saturation curve
# ----------------
Ts = np.linspace(Tmin, Tc, 1000)
ps = CP.CoolProp.PropsSI('P','T',Ts,'Q',0,fluidname)

# ----------------
# Contour of Z and Gamma
# ----------------
n = 200 # number of points
x = np.linspace(Tmin, Tmax,n)
y = np.linspace(1e5, 1e7,n)
X,Y = np.meshgrid(x,y)
Z =  CP.CoolProp.PropsSI('Z','T',X,'P',Y,fluidname)
Gamma =  CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','T',X,'P',Y,fluidname)
#print(Z.shape)
levels = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
cp = plt.contour(X, Y, Z, levels, colors='black', linestyles='dashed')
plt.clabel(cp, inline=True,  fontsize=10)
plt.contourf(X, Y, Gamma, [0.4,0.6,0.8,1.0,1.2,1.4,1.6], cmap='rainbow')
plt.colorbar()
# ------
# Labels
# ------

plt.plot(Ts,ps,'k',lw = lw, solid_capstyle = 'round', label = "LVS")
ax.legend(loc=4) # 2 means left top
# Critical lines
plt.axvline(Tc, dashes = [2, 2])
plt.axhline(pc, dashes = [2, 2])

"""
test points
"""

nc = 10
colors = plt.cm.tab20(np.linspace(0, 1, nc))

z9_p = [1.55e6,1.30E6,  1.05E6, 8.01E5, 5.51E5, 3.01E5]
z9_t = [673, 654.39, 624.85, 585.46, 531.30, 449.23]
plt.plot(z9_t,z9_p,'o' ,color=colors[0], lw = lw)

z8_p = [2.52e6,2.12E6, 1.72E6, 1.32E6, 9.20E5, 5.21E5]
z8_t = [643.45, 623.76, 598.05, 564.13, 515.98, 443.76]
plt.plot(z8_t,z8_p,'o' ,color=colors[1],lw = lw)

z7_p = [2.91E6, 2.51E6, 2.11E6, 1.71E6, 1.31E6, 9.08E5 ]
z7_t = [610.08, 595.31, 576.16, 551.54, 519.26, 474.40]
plt.plot(z7_t,z7_p,'o' , color=colors[2],lw = lw)

z6_p = [5e6, 4.3e6,  3.6e6, 2.9e6, 2.2e6, 1.5e6]
z6_t = [619.38, 612.27, 599.69, 581.08, 553.733, 511.60]
plt.plot(z6_t,z6_p,'o' , color=colors[3],lw = lw)

z5_p = [5e6, 4.35e6, 3.70e6, 3.05E6, 2.40E6, 1.75E6]
z5_t = [589.29, 588.20, 580.54, 567.41, 546.62, 515.98]
plt.plot(z5_t,z5_p,'o' ,color=colors[4],lw = lw)

z4_p = [4e6,   3.6e6, 3.2e6, 2.8e6, 2.4e6, 2e6]
z4_t = [563.03,  561.94, 557.01, 548.81, 539.86, 523.09]
plt.plot(z4_t,z4_p,'o' , color=colors[5],lw = lw)

plt.ylim(1e5,1e7)
plt.gca().set_yscale('log')
plt.gca().set_xlim(Tmin, Tmax)
plt.ylabel('Pressure [Pa]')
plt.xlabel('Temperature [K]')
plt.title('Contour of Z and $\Gamma$ for siloxane MM')
plt.tight_layout()
fig.savefig("files/mm_z_Contour_PT.eps")
print("plotcontour.py called")