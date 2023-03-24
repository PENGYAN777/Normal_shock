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
from newIOpairs import TGfromZP,PGfromZT


"""
0. fluid property
"""
# print("--------------------Fluid info ------------------")
fluidname = "MM"
print("Fluid name:", fluidname)
R = CP.CoolProp.PropsSI("gas_constant",fluidname)
W = CP.CoolProp.PropsSI("molar_mass",fluidname)
Rs = R/W
print("spefific ags constant: J/Kg/K", Rs)
Tc =  CP.CoolProp.PropsSI("Tcrit",fluidname)
print("critical temperature[K]:", Tc)
Pc =  CP.CoolProp.PropsSI("pcrit",fluidname)
print("critical pressure[Pa]:", Pc)
dc = CP.CoolProp.PropsSI('Dmass','P',Pc,'T',Tc,fluidname) 

Zt = np.arange(0.70, 0.95, 0.02).tolist()
r_m = np.zeros(len(Zt)) # M2/M1
r_p = np.zeros(len(Zt)) # P2/P1
r_t = np.zeros(len(Zt)) # T2/T1
r_d = np.zeros(len(Zt)) # D2/D1
r_pt = np.zeros(len(Zt)) # PT2/PT1
r_y = np.zeros(len(Zt)) # shock losee

for j in range(0,len(Zt),1):
    """
    1. total conditions
    """
    
    print("--------------------Total conditions ------------------")
    tt = 550 # total pressure
    zt = Zt[j]
    pt, gt = PGfromZT(zt,tt)
    s = CP.CoolProp.PropsSI('Smass','P',pt,'T',tt,fluidname) 
    ht = CP.CoolProp.PropsSI('Hmass','P',pt,'T',tt,fluidname) 
    
    """
    1.1. compute isentropic relationship
    """
    
    p = np.linspace(pt*0.1,pt,1000) #
    p = pd.Series(p)
    h = np.zeros(p.size) # enthalpy
    u = np.zeros(p.size) # velocity
    c = np.zeros(p.size) # sound speed
    m = np.zeros(p.size) # Mach number
    d = np.zeros(p.size) # density
    t = np.zeros(p.size) # temperature
    for i in p.index:
        if abs(p[i]-Pc)<0.01*Pc:
            p[i] = 0.99*Pc
        d[i] = CP.CoolProp.PropsSI('Dmass','P',p[i],'Smass',s,fluidname) 
        t[i] = CP.CoolProp.PropsSI('T','P',p[i],'Smass',s,fluidname) 
        h[i] = CP.CoolProp.PropsSI('Hmass','P',p[i],'Smass',s,fluidname) 
        # h[i] = CP.CoolProp.PropsSI('Hmass','T',t[i],'P',p[i],fluidname)
        u[i] = math.sqrt(abs(2*(ht-h[i])))
        # c[i] = CP.CoolProp.PropsSI('A','P',p[i],'T',t[i],fluidname) 
        c[i] = CP.CoolProp.PropsSI('A','P',p[i],'Smass',s,fluidname) 
        m[i] = u[i]/c[i]
    
    """
    1.2. find pre-shock condition, M1
    """
    # print("--------------------Pre-shock states ------------------")
    M1 = 1.5
    # print('M1: ',M1)
    # print("index for required condition:",np.argmin(abs(m-M1)))
    d1 = d[np.argmin(abs(m-M1))]
    T1 = t[np.argmin(abs(m-M1))]
    P1 = p[np.argmin(abs(m-M1))]
    # c1 = c[np.argmin(abs(m-M1))]
    u1 = u[np.argmin(abs(m-M1))]
    # h1 = h[np.argmin(abs(m-M1))]
    # print('P1[Pa]: ',P1)
    # print('T1[Pa]: ',T1)
    
    
    """
    2. input pre-shock conditions
    """ 
    
    """
    3. compute post-shock properites
    """
    # print("--------------------Post-shock states ------------------")
    p2 = np.linspace(P1*1.1,P1*3 ,1000) # post-shock pressure 
    p2 = pd.Series(p2)
    u2 = np.zeros(p2.size) 
    d2 = np.zeros(p2.size) 
    diff = np.zeros(p2.size) + 1e3
     
    for i in p2.index:
        if abs(p2[i]-Pc)<0.01*Pc:
            p2[i] = 0.99*Pc
        u2[i] = (P1+d1*u1*u1-p2[i])/d1/u1
        if u2[i]>0:
            d2[i] =  d1*u1/u2[i]
            diff[i] = ht - 0.5*u2[i]*u2[i] - CP.CoolProp.PropsSI('Hmass','P',p2[i],'Dmass',d2[i],fluidname) 
            if abs(diff[i])<1e-4:
                break
    # print("min index:", np.argmin(abs(diff)))
    P2 = p2[np.argmin(abs(diff))]
    D2 = d2[np.argmin(abs(diff))]    
    T2 = CP.CoolProp.PropsSI('T','P',P2,'Dmass',D2,fluidname) 
    c2 = CP.CoolProp.PropsSI('A','P',P2,'Dmass',D2,fluidname) 
    U2 = u2[np.argmin(abs(diff))]    
    M2 = U2/c2
    h2 = CP.CoolProp.PropsSI('Hmass','P',P2,'Dmass',D2,fluidname)
    ht2 = h2 + 0.5*U2*U2
    # print("ht2:", htotal2)
    print("(ht2-ht1)/ht1(%):", (ht2-ht)/ht*100)
    # print("M2:", M2)
    # print("P2/P1:", P2/P1)
    # print("T2/T1:", T2/T1)
    # print("D2/D1:", D2/d1)
    s2 = CP.CoolProp.PropsSI('Smass','P',P2,'Dmass',D2,fluidname) 
    Pt2 = CP.CoolProp.PropsSI('P','Smass',s2,'Hmass',ht2,fluidname) 
    r_m[j] = M2/M1
    r_p[j] = P2/P1
    r_t[j] = T2/T1
    r_d[j] = D2/d1
    r_pt[j] = Pt2/pt
    r_y[j] = (pt-Pt2)/(pt-P1)*100
    # Tt2 = CP.CoolProp.PropsSI('T','Smass',s2,'Hmass',ht2,fluidname) 
    # Y = (pt-Pt2)/(pt-P1)*100

"""
4. write into csv file
"""   

pd.DataFrame(Zt).to_csv('data.csv', index_label = "Index", header  = ['Zt']) 
data = pd.read_csv("data.csv", ",")
# append new columns
D =pd.DataFrame({'P2/P1': r_p, 'T2/T1': r_t, 'D2/D1': r_d,'M2/M1': r_m, 'Pt2/Pt1': r_pt,'Y': r_y,})
newData = pd.concat([data, D], join = 'outer', axis = 1)
# save newData in csv file
# newData.to_csv("m4sh.csv")
newData.to_csv("data.csv")
