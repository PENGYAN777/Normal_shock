#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:35:29 2022

@author: yan

compute normal shock relatiomnships for CoolProp EOS
input: upstream P,T,M.

"""

from scipy.optimize import fsolve
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import CoolProp as CP


"""
0. fluid property
"""
print("--------------------Fluid info ------------------")
fluidname = "MM"
print("Fluid name:", fluidname)
R = CP.CoolProp.PropsSI("gas_constant",fluidname)
# print("universal gas constant:  J/mol/K", R)
W = CP.CoolProp.PropsSI("molar_mass",fluidname)
# print("molar mass: kg/mol", W)
Rs = R/W
print("spefific ags constant: J/Kg/K", Rs)
Tc =  CP.CoolProp.PropsSI("Tcrit",fluidname)
print("critical temperature[K]:", Tc)
Pc =  CP.CoolProp.PropsSI("pcrit",fluidname)
print("critical pressure[Pa]:", Pc)
dc = CP.CoolProp.PropsSI('Dmass','P',Pc,'T',Tc,fluidname) 

"""
1. total conditions
"""

print("--------------------Total conditions ------------------")
pt = 1e6 # total pressure
tt = 225+273.15
zt = CP.CoolProp.PropsSI('Z','P',pt,'T',tt,fluidname) 
print('Pt[Pa]: ',pt)
print('Tt[K]: ',tt)
print('Zt: ',zt)
# dt = CP.CoolProp.PropsSI('Dmass','P',pt,'T',tt,fluidname) 
s = CP.CoolProp.PropsSI('Smass','P',pt,'T',tt,fluidname) 
ht = CP.CoolProp.PropsSI('Hmass','P',pt,'T',tt,fluidname) 
"""
1.1. compute isentropic relationship
"""

p = np.linspace(pt*0.1,pt,1000) #

p = pd.Series(p)
Z = np.zeros(p.size) # P/rho RT
h = np.zeros(p.size) # enthalpy
u = np.zeros(p.size) # velocity
v = np.zeros(p.size) # specific volume
c = np.zeros(p.size) # sound speed
m = np.zeros(p.size) # Mach number
d = np.zeros(p.size) # density
t = np.zeros(p.size) # temperature
# g = np.zeros(p.size) # Gamma

for i in p.index:
    if abs(p[i]-Pc)<0.01*Pc:
        p[i] = 0.99*Pc
    d[i] = CP.CoolProp.PropsSI('Dmass','P',p[i],'Smass',s,fluidname) 
    t[i] = CP.CoolProp.PropsSI('T','P',p[i],'Smass',s,fluidname) 
    Z[i] = CP.CoolProp.PropsSI('Z','P',p[i],'T',t[i],fluidname) 
    # g[i] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','P',p[i],'Smass',s,fluidname) 
    h[i] = CP.CoolProp.PropsSI('Hmass','T',t[i],'P',p[i],fluidname)
    u[i] = math.sqrt(abs(2*(ht-h[i])))
    c[i] = CP.CoolProp.PropsSI('A','P',p[i],'T',t[i],fluidname) 
    m[i] = u[i]/c[i]

"""
1.2. find pre-shock condition, M1
"""
print("--------------------Pre-shock states ------------------")
M1 = 1.5
print('M1: ',M1)
print("index for required condition:",np.argmin(abs(m-M1)))
d1 = d[np.argmin(abs(m-M1))]
T1 = t[np.argmin(abs(m-M1))]
P1 = p[np.argmin(abs(m-M1))]
c1 = c[np.argmin(abs(m-M1))]
u1 = u[np.argmin(abs(m-M1))]
h1 = h[np.argmin(abs(m-M1))]
print('P1[Pa]: ',P1)
print('T1[Pa]: ',T1)


"""
2. input pre-shock conditions
"""
# d1 = CP.CoolProp.PropsSI('Dmass','P',P1,'T',T1,fluidname) # pre-shock density
Z1 = CP.CoolProp.PropsSI('Z','P',P1,'T',T1,fluidname) 
G1 = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','P',P1,'T',T1,fluidname) 
print("Z1,Gamma1:", Z1,G1)
s1 = CP.CoolProp.PropsSI('Smass','P',P1,'T',T1,fluidname) 
print("st:", s)
print("s1:", s1)
cv = CP.CoolProp.PropsSI('Cvmass','P',P1,'T',T1,fluidname) 
cp = CP.CoolProp.PropsSI('Cpmass','P',P1,'T',T1,fluidname) 
print("gamma:", cp/cv)
ht1 = h1 + 0.5*u1*u1
print("ht:", ht)
print("h1:", ht1)
# compute total pressure/temperature
Pt1 = CP.CoolProp.PropsSI('P','Smass',s1,'Hmass',ht1,fluidname) 
Tt1 = CP.CoolProp.PropsSI('T','Smass',s1,'Hmass',ht1,fluidname) 
htotal1 = CP.CoolProp.PropsSI('Hmass','P',Pt1,'T',Tt1,fluidname) 
print("ht1:", htotal1)

"""
3. compute post-shock properites
"""
print("--------------------Post-shock states ------------------")
p2 = np.linspace(P1*1.1,P1*3.0 ,1000) # post-shock pressure 
p2 = pd.Series(p2)
u2 = np.zeros(p2.size) 
d2 = np.zeros(p2.size) 
diff = np.zeros(p2.size) 
 
for i in p2.index:
    if abs(p2[i]-Pc)<0.01*Pc:
        p2[i] = 0.99*Pc
    u2[i] = (P1+d1*u1*u1-p2[i])/d1/u1
    if u2[i]>0:
        d2[i] =  d1*u1/u2[i]
        diff[i] = ht1 - 0.5*u2[i]*u2[i] - CP.CoolProp.PropsSI('Hmass','P',p2[i],'Dmass',d2[i],fluidname) 
print("min index:", np.argmin(abs(diff)))
P2 = p2[np.argmin(abs(diff))]
D2 = d2[np.argmin(abs(diff))]    
T2 = CP.CoolProp.PropsSI('T','P',P2,'Dmass',D2,fluidname) 
c2 = CP.CoolProp.PropsSI('A','P',P2,'T',T2,fluidname) 
U2 = u2[np.argmin(abs(diff))]    
M2 = U2/c2
h2 = CP.CoolProp.PropsSI('Hmass','P',P2,'Dmass',D2,fluidname) 
htotal2 = h2 + 0.5*U2*U2
print("ht2:", htotal2)
print("(ht2-ht1)/ht1(%):", (htotal2-htotal1)/htotal1*100)
print("M2:", M2)
print("P2/P1:", P2/P1)
print("T2/T1:", T2/T1)
print("D2/D1:", D2/d1)
s2 = CP.CoolProp.PropsSI('Smass','P',P2,'T',T2,fluidname) 
Pt2 = CP.CoolProp.PropsSI('P','Smass',s2,'Hmass',ht1,fluidname) 
Tt2 = CP.CoolProp.PropsSI('T','Smass',s2,'Hmass',ht1,fluidname) 
print("Pt2/Pt1:", Pt2/Pt1)
print("Y(%) = (Pt1-Pt2)/(Pt1-P1)*100:", (Pt1-Pt2)/(Pt1-P1)*100)
