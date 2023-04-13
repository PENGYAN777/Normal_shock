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
from newIOpairs import TGfromZP,PGfromZT,ZPfromTG


"""
0. fluid property
"""
# print("--------------------Fluid info ------------------")
fluidname = "D4"
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

G1 = np.arange(0.6, 0.92, 0.05).tolist()
M1 = 1.5
for k in range(0,len(G1),1):
    print("---------------G1------------:", G1[k])
    T1_list =np.linspace(700-200*G1[k], 50*G1[k]+570 , 6)
    T1_list = T1_list.tolist()
    r_m = np.zeros(len(T1_list)) # M2/M1
    r_p = np.zeros(len(T1_list)) # P2/P1
    r_t = np.zeros(len(T1_list)) # T2/T1
    r_d = np.zeros(len(T1_list)) # D2/D1
    r_pt = np.zeros(len(T1_list)) # PT2/PT1
    r_y = np.zeros(len(T1_list)) # shock losee
    Diff = np.zeros(len(T1_list)) # shock losee
    pfs= np.zeros(len(T1_list)) # freestram pressure
    tfs= np.zeros(len(T1_list)) # freestram temperature
    G_value = np.zeros(len(T1_list)) # 
    for j in range(0,len(T1_list),1):
        """
        1. pre-shock conditions
        """
        
        print("--------------------pre-shock conditions ------------------")
        T1 = T1_list[j]
        Z1, P1 = ZPfromTG(T1,G1[k])
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
        p2 = np.linspace(P1*1.5,P1*(2.18 + 0.2*k/1.5) ,1000) # post-shock pressure 
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
                    print('diff is:', diff[i])
                    break
        # print("min index:", np.argmin(abs(diff)))
        P2 = p2[np.argmin(abs(diff))]
        D2 = d2[np.argmin(abs(diff))]    
        T2 = CP.CoolProp.PropsSI('T','P',P2,'Dmass',D2,fluidname) 
        c2 = CP.CoolProp.PropsSI('A','P',P2,'T|gas',T2,fluidname) 
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
        G_value[j] = G1[k]
        # Tt2 = CP.CoolProp.PropsSI('T','Smass',s2,'Hmass',ht2,fluidname) 
        # Y = (pt-Pt2)/(pt-P1)*100
    
    """
    4. write into csv file
    """   
    data = 'data' + str(k) + '.csv'
    pd.DataFrame(pfs).to_csv(data, index_label = "Index", header  = ['P1']) 
    result = pd.read_csv(data, ",")
    # append new columns
    D =pd.DataFrame({'P2/P1': r_p, 'T2/T1': r_t, 'D2/D1': r_d,'M2/M1': r_m, 'Pt2/Pt1': r_pt,'Y': r_y, 
                     'diff': Diff, 'T1': tfs, 'G1': G_value, })
    newData = pd.concat([result, D], join = 'outer', axis = 1)
    # save newData in csv file
    # newData.to_csv("m4sh.csv")
    newData.to_csv(data)