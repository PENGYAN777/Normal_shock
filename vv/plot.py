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
Yc= pd.read_csv("Yc.csv", ",", header=None)


# n = 10
# colors = plt.cm.tab20(np.linspace(0, 1, n))
fig1, ax = plt.subplots(3)
fig1.suptitle('Comparison between theoretical solutions and previous work')
# fig1 = plt.figure( dpi=300)
lwh = 2
# ax[0] = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax[0].plot(data.iloc[:,2]/1e5 , data.iloc[:,-2] , 'ko', lw=lwh, label="Experiment data (Conti et.al 2022)")
ax[0].plot(data.iloc[:,2]/1e5 , Yc.iloc[:,-1] , 'bo', lw=lwh, label="Calculated results (Conti et.al 2022)")
ax[0].plot(data.iloc[:,2]/1e5 , data.iloc[:,-1] , 'ro', lw=lwh, label="Theoretical solutions")


ax[0].set_xlabel('$P_t$ $[bar]$',fontsize=12)
ax[0].set_ylabel('Y(%)',fontsize=12) 
# axes.set_title('Mach number vs $\\rho/\\rho_t$',fontsize=14)
ax[0].set_ylim([15, 25])
ax[0].legend(loc=0 , prop={'size': 7}) # 

ax[1].plot(data.iloc[:,2]/1e5 , data.iloc[:,4] , 'ko', lw=lwh)
ax[1].set_xlabel('$P_t$ $[bar]$',fontsize=12)
ax[1].set_ylabel('$M_1$',fontsize=12) 
# axes.set_title('Mach number vs $\\rho/\\rho_t$',fontsize=14)
ax[1].set_ylim([1.4, 1.6])

# ax[2] = fig1.add_axes([0.15, 0.15, 0.2, 0.2]) #size of figure
ax[2].plot(data.iloc[:,2]/1e5 , data.iloc[:,3]-273.15 , 'ko')
ax[2].set_xlabel('$P_t$ $[bar]$',fontsize=12)
ax[2].set_ylabel('$T_t$ [K]',fontsize=12) 
# axes.set_title('Mach number vs $\\rho/\\rho_t$',fontsize=14)
ax[2].set_ylim([200, 240])
# ax[1].legend(loc=0 , prop={'size': 10}) # 

fig1.savefig("vv.eps")



# fig3 = plt.figure( dpi=300)
# lwh = 2
# axes = fig3.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(data.iloc[:,2]/1e5 , data.iloc[:,4] , 'ko', lw=lwh, label="Experiment data")
# axes.set_xlabel('$P_t$ $[bar]$',fontsize=12)
# axes.set_ylabel('$M_1$',fontsize=12) 
# # axes.set_title('Mach number vs $\\rho/\\rho_t$',fontsize=14)
# axes.set_ylim([1.35, 1.65])
# axes.legend(loc=0 , prop={'size': 10}) # 