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
T1_list = np.arange(540, 580, 10).tolist()
nc = len(T1_list)
colors = plt.cm.tab20(np.linspace(0, 1, nc))

"""
plot 
"""
nn = 10
fig1 = plt.figure( dpi=300)
lw = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax2 = axes.twinx()
pmax = np.zeros(nn)
pmin = np.zeros(nn)+100
for k in range(0,len(T1_list),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)
    lb = 'T1=' + str(T1_list[k]) + '(K)'
    axes.plot(z.iloc[:,2],z.iloc[:,3],'o', color=colors[k],  lw = lw, label = lb)
    for i in range(10):
        pmax[i] = max(pmax[i],z.iloc[i,3])
        pmin[i] = min(pmin[i],z.iloc[i,3])
    
ax2.plot(z.iloc[:,2] , (pmax-pmin)/pmax*100 , 'k*', lw=lw)
ax2.set_ylabel('diff(%)',fontsize=12) 
ax2.set_ylim([0, 3])
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.set_xlim([0.5, 1.0])
axes.set_ylim([1.6, 2.3])
axes.set_xlabel('$Z_1$',fontsize=12)
axes.set_ylabel('$P_2/P_1$',fontsize=12) 
axes.legend(loc=8 , prop={'size': 10}) # 
axes.set_title('$P_2/P_1$ vs $Z_t$',fontsize=14)
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig1.savefig("files/nicfd_mm_rp.eps")

fig2 = plt.figure( dpi=300)
lw = 2
axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax2 = axes.twinx()
tmax = np.zeros(nn)
tmin = np.zeros(nn)+100
for k in range(0,len(T1_list),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)
    lb = 'T1=' + str(T1_list[k]) + '(K)'
    axes.plot(z.iloc[:,2],z.iloc[:,4],'o', color=colors[k],  lw = lw, label = lb)
    for i in range(10):
        tmax[i] = max(tmax[i],z.iloc[i,4])
        tmin[i] = min(tmin[i],z.iloc[i,4])
    
ax2.plot(z.iloc[:,2] , (tmax-tmin)/tmax*100 , 'k*', lw=lw)
ax2.set_ylabel('diff(%)',fontsize=12) 
# ax2.set_ylim([0, 3])
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.set_xlim([0.5, 1.0])
axes.set_xlabel('$Z_1$',fontsize=12)
axes.set_ylabel('$T_2/T_1$',fontsize=12) 
axes.legend(loc=0 , prop={'size': 10}) # 
axes.set_title('$T_2/T_1$ vs $Z_t$',fontsize=14)
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig2.savefig("files/nicfd_mm_rt.eps")

fig3 = plt.figure( dpi=300)
lw = 2
axes = fig3.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax2 = axes.twinx()
dmax = np.zeros(nn)
dmin = np.zeros(nn)+100
for k in range(0,len(T1_list),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)
    lb = 'T1=' + str(T1_list[k]) + '(K)'
    axes.plot(z.iloc[:,2],z.iloc[:,5],'o', color=colors[k],  lw = lw, label = lb)
    for i in range(10):
        dmax[i] = max(dmax[i],z.iloc[i,5])
        dmin[i] = min(dmin[i],z.iloc[i,5])
    
ax2.plot(z.iloc[:,2] , (dmax-dmin)/dmax*100 , 'k*', lw=lw)
ax2.set_ylabel('diff(%)',fontsize=12) 
# ax2.set_ylim([0, 3])
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.set_xlim([0.5, 1.0])
axes.set_xlabel('$Z_1$',fontsize=12)
axes.set_ylabel('$\\rho_2/\\rho_1$',fontsize=12) 
axes.legend(loc=8 , prop={'size': 10}) # 
axes.set_title('$\\rho_2/\\rho_1$ vs $Z_t$',fontsize=14)
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig3.savefig("files/nicfd_mm_rd.eps")

fig4 = plt.figure( dpi=300)
lw = 2
axes = fig4.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax2 = axes.twinx()
mmax = np.zeros(nn)
mmin = np.zeros(nn)+100
for k in range(0,len(T1_list),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)
    lb = 'T1=' + str(T1_list[k]) + '(K)'
    axes.plot(z.iloc[:,2],z.iloc[:,6],'o', color=colors[k],  lw = lw, label = lb)
    for i in range(10):
        mmax[i] = max(mmax[i],z.iloc[i,6])
        mmin[i] = min(mmin[i],z.iloc[i,6])
    
ax2.plot(z.iloc[:,2] , (mmax-mmin)/mmax*100 , 'k*', lw=lw)
ax2.set_ylabel('diff(%)',fontsize=12) 
# ax2.set_ylim([0, 3])
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.set_xlim([0.5, 1.0])
axes.set_xlabel('$Z_1$',fontsize=12)
axes.set_ylabel('$M_2/M_1$',fontsize=12) 
axes.legend(loc=5 , prop={'size': 10}) # 
axes.set_title('$M_2/M_1$ vs $Z_t$',fontsize=14)
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig4.savefig("files/nicfd_mm_rm.eps")

fig5 = plt.figure( dpi=300)
lw = 2
axes = fig5.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax2 = axes.twinx()
ptmax = np.zeros(nn)
ptmin = np.zeros(nn)+100
for k in range(0,len(T1_list),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)
    lb = 'T1=' + str(T1_list[k]) + '(K)'
    axes.plot(z.iloc[:,2],z.iloc[:,7],'o', color=colors[k],  lw = lw, label = lb)
    for i in range(10):
        ptmax[i] = max(ptmax[i],z.iloc[i,7])
        ptmin[i] = min(ptmin[i],z.iloc[i,7])
    
ax2.plot(z.iloc[:,2] , (ptmax-ptmin)/ptmax*100 , 'k*', lw=lw)
ax2.set_ylabel('diff(%)',fontsize=12) 
ax2.set_ylim([0, 5])
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.set_xlim([0.5, 1.0])
axes.set_xlabel('$Z_1$',fontsize=12)
axes.set_ylabel('$P_{t2}/P_{t1}$',fontsize=12) 
axes.legend(loc=0 , prop={'size': 10}) # 
axes.set_title('$P_{t2}/P_{t1}$ vs $Z_t$',fontsize=14)
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig5.savefig("files/nicfd_mm_rpt.eps")

fig6 = plt.figure( dpi=300)
lw = 2
axes = fig6.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax2 = axes.twinx()
ymax = np.zeros(nn)
ymin = np.zeros(nn)+100
for k in range(0,len(T1_list),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)
    lb = 'T1=' + str(T1_list[k]) + '(K)'
    axes.plot(z.iloc[:,2],z.iloc[:,8],'o', color=colors[k],  lw = lw, label = lb)
    for i in range(10):
        ymax[i] = max(ymax[i],z.iloc[i,8])
        ymin[i] = min(ymin[i],z.iloc[i,8])
    
ax2.plot(z.iloc[:,2] , (ymax-ymin)/ymax*100 , 'k*', lw=lw)
ax2.set_ylabel('diff(%)',fontsize=12) 
# ax2.set_ylim([0, 3])
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.set_xlim([0.5, 1.0])
axes.set_xlabel('$Z_1$',fontsize=12)
axes.set_ylabel('Y(%)',fontsize=12) 
axes.legend(loc=0 , prop={'size': 10}) # 
axes.set_title('Y(%) vs $Z_t$',fontsize=14)
axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig6.savefig("files/nicfd_mm_ry.eps")





# fig10 = plt.figure( dpi=300)
# lw = 2
# axes = fig10.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# for k in range(0,len(T1_list),1):
#     data = 'data' + str(k) + '.csv'
#     z = pd.read_csv(data, ",", skiprows=0)
#     lb = 'T1=' + str(T1_list[k]) + '(K)'
#     axes.plot(z.iloc[:,2],z.iloc[:,-3],'o', color=colors[k],  lw = lw, label = lb)
# axes.set_xlabel('$Z_1$',fontsize=12)
# axes.set_ylabel('$(h_{t2}-h_{t1})/h_{t2}(\%)$',fontsize=12) 
# axes.legend(loc=0 , prop={'size': 10}) # 
# axes.set_title('$(h_{t2}-h_{t1})/h_{t1}(\%)$ vs $Z_t$',fontsize=14)
# axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# # axes.set_ylim([0, 0.1])
# # fig1.savefig("files/nicfd_mm_rh.eps")


