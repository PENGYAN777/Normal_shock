#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 16:36:16 2023

@author: yan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

data= pd.read_csv("result.csv", ",", skiprows=0)


# n = 10
# colors = plt.cm.tab20(np.linspace(0, 1, n))
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(data.iloc[:,2]/1e5 , data.iloc[:,-2]*1.05 , 'k', lw=lwh, label="$\pm$ 5% bar")
axes.plot(data.iloc[:,2]/1e5 , data.iloc[:,-2]*0.95 , 'k', lw=lwh)
axes.plot(data.iloc[:,2]/1e5 , data.iloc[:,-2] , 'ko', lw=lwh, label="Experiment data")
axes.plot(data.iloc[:,2]/1e5 , data.iloc[:,-1] , 'ro', lw=lwh, label="theoretical solutions")


axes.set_xlabel('$P_t$ $[bar]$',fontsize=12)
axes.set_ylabel('Y(%)',fontsize=12) 
# axes.set_title('Mach number vs $\\rho/\\rho_t$',fontsize=14)
axes.set_ylim([10, 20])
axes.legend(loc=0 , prop={'size': 10}) # 
fig1.savefig("vv.eps")

fig2 = plt.figure( dpi=300)
lwh = 2
axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(data.iloc[:,2]/1e5 , data.iloc[:,3]-273.15 , 'ko', lw=lwh, label="Experiment data")
axes.set_xlabel('$P_t$ $[bar]$',fontsize=12)
axes.set_ylabel('$T_t$ [K]',fontsize=12) 
# axes.set_title('Mach number vs $\\rho/\\rho_t$',fontsize=14)
axes.set_ylim([160, 240])
axes.legend(loc=0 , prop={'size': 10}) # 

fig3 = plt.figure( dpi=300)
lwh = 2
axes = fig3.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(data.iloc[:,2]/1e5 , data.iloc[:,4] , 'ko', lw=lwh, label="Experiment data")
axes.set_xlabel('$P_t$ $[bar]$',fontsize=12)
axes.set_ylabel('Mach',fontsize=12) 
# axes.set_title('Mach number vs $\\rho/\\rho_t$',fontsize=14)
axes.set_ylim([1.35, 1.65])
axes.legend(loc=0 , prop={'size': 10}) # 