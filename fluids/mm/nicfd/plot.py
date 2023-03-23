#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 14:59:53 2023

@author: yan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

"""
read csv file
"""

data = pd.read_csv("data.csv", ",", skiprows=0)

"""
plot 
"""
n = 10
colors = plt.cm.tab20(np.linspace(0, 1, n))

fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(data.iloc[:,2] , data.iloc[:,3] , 'ko', lw=lwh, label="$P_2/P_1$")

axes.set_xlabel('$Z_t$',fontsize=12)
axes.set_ylabel('$P_2/P_1$',fontsize=12) 
axes.set_title('$P_2/P_1$ vs $Z_t$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig1.savefig("files/nicfd_mm_rp.eps")

fig2 = plt.figure( dpi=300)
lwh = 2
axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(data.iloc[:,2] , data.iloc[:,4] , 'ko', lw=lwh, label="$T_2/T_1$")

axes.set_xlabel('$Z_t$',fontsize=12)
axes.set_ylabel('$T_2/T_1$',fontsize=12) 
axes.set_title('$T_2/T_1$ vs $Z_t$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig2.savefig("files/nicfd_mm_rT.eps")

fig3 = plt.figure( dpi=300)
lwh = 2
axes = fig3.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(data.iloc[:,2] , data.iloc[:,5] , 'ko', lw=lwh, label="$\\rho_2/\\rho_1$")

axes.set_xlabel('$Z_t$',fontsize=12)
axes.set_ylabel('$\\rho_2/\\rho_1$',fontsize=12) 
axes.set_title('$\\rho_2/\\rho_1$ vs $Z_t$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig3.savefig("files/nicfd_mm_rd.eps")

fig4 = plt.figure( dpi=300)
lwh = 2
axes = fig4.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(data.iloc[:,2] , data.iloc[:,6] , 'ko', lw=lwh, label="$M_2/M_1$")

axes.set_xlabel('$Z_t$',fontsize=12)
axes.set_ylabel('$M_2/M_1$',fontsize=12) 
axes.set_title('$M_2/M_1$ vs $Z_t$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig4.savefig("files/nicfd_mm_rm.eps")

fig5 = plt.figure( dpi=300)
lwh = 2
axes = fig5.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(data.iloc[:,2] , data.iloc[:,7] , 'ko', lw=lwh, label="$Pt_2/Pt_1$")

axes.set_xlabel('$Z_t$',fontsize=12)
axes.set_ylabel('$Pt_2/Pt_1$',fontsize=12) 
axes.set_title('$Pt_2/Pt_1$ vs $Z_t$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig5.savefig("files/nicfd_mm_rpt.eps")

fig6 = plt.figure( dpi=300)
lwh = 2
axes = fig6.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(data.iloc[:,2] , data.iloc[:,8] , 'ko', lw=lwh, label="Y(%)")

axes.set_xlabel('$Z_t$',fontsize=12)
axes.set_ylabel('Y(%)',fontsize=12) 
axes.set_title('Y(%)',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig5.savefig("files/nicfd_mm_ry.eps")