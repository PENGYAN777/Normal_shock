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

Z1 = np.arange(0.50, 0.82, 0.05).tolist()
r_m = np.zeros(len(Z1)) # M2/M1
r_p = np.zeros(len(Z1)) # P2/P1
r_t = np.zeros(len(Z1)) # T2/T1
r_d = np.zeros(len(Z1)) # D2/D1
r_pt = np.zeros(len(Z1)) # PT2/PT1
r_y = np.zeros(len(Z1)) # shock loss
Diff = np.zeros(len(Z1)) # tolerance
pfs= np.zeros(len(Z1)) # freestram pressure
tfs= np.zeros(len(Z1)) # freestram temperature

P1_list = np.arange(1.1*Pc, 1.5*Pc, 0.1*Pc).tolist()
M1 = 1.5
for k in range(0,len(P1_list),1):
    print("--------------------P1/Pc------------------",P1_list[k]/Pc)
    for j in range(0,len(Z1),1):
        """
        1. pre-shock conditions
        """
        P1 = P1_list[k] # total pressure
        z1 = Z1[j]
        T1, g1 = TGfromZP(z1,P1)
        d1 = CP.CoolProp.PropsSI('Dmass','P',P1,'T',T1,fluidname) 
        s = CP.CoolProp.PropsSI('Smass','P',P1,'T',T1,fluidname) 
        h1 = CP.CoolProp.PropsSI('Hmass','P',P1,'T',T1,fluidname) 
        c1 = CP.CoolProp.PropsSI('A','P',P1,'T',T1,fluidname) 
        u1 = c1*M1
        ht = h1 + 0.5*u1*u1
        Pt1 = CP.CoolProp.PropsSI('P','Smass',s,'Hmass',ht,fluidname) 
        """
        2. input pre-shock conditions
        """ 
        
        """
        3. compute post-shock properites
        """
        # print("--------------------Post-shock states ------------------")
        p2 = np.linspace(P1*1.1,P1*(1.8 + 0.1*j), 1000) # post-shock pressure 
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
                if abs(diff[i])<1:
                    print('diff is: ', diff[i])
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
        # print("(ht2-ht1)/ht1(%):", (ht2-ht)/ht*100)
        Diff[j] = (ht2-ht)/ht*100
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
        r_pt[j] = Pt2/Pt1
        r_y[j] = (Pt1-Pt2)/(Pt1-P1)*100
        pfs[j] = P1
        tfs[j] = T1
        # Tt2 = CP.CoolProp.PropsSI('T','Smass',s2,'Hmass',ht2,fluidname) 
        # Y = (pt-Pt2)/(pt-P1)*100
    
    """
    4. write into csv file
    """   
    data = 'data' + str(k) + '.csv'
    pd.DataFrame(Z1).to_csv(data, index_label = "Index", header  = ['Z1']) 
    result = pd.read_csv(data, ",")
    # append new columns
    D =pd.DataFrame({'P2/P1': r_p, 'T2/T1': r_t, 'D2/D1': r_d,'M2/M1': r_m, 'Pt2/Pt1': r_pt,'Y': r_y, 
                     'diff': Diff, 'P1': pfs, 'T1': tfs, })
    newData = pd.concat([result, D], join = 'outer', axis = 1)
    # save newData in csv file
    # newData.to_csv("m4sh.csv")
    newData.to_csv(data)
