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
G1 = np.arange(0.6, 0.92, 0.05).tolist()
# for k in range(0,len(Z1),1):
#     data = 'data' + str(k) + '.csv'
#     z = pd.read_csv(data, ",", skiprows=0)

    

"""
plot 
"""
nc = len(G1)
colors = plt.cm.tab20(np.linspace(0, 1, nc))

fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax2 = axes.twinx()
for k in range(0,len(G1),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)
    axes.plot(z.iloc[:,-1] , z.iloc[:,3] , 'o',color=colors[k], lw=lwh)
    max_value = z.iloc[:,3][np.argmax(abs(z.iloc[:,3]))]
    min_value = z.iloc[:,3][np.argmin(abs(z.iloc[:,3]))]
    diff = (max_value-min_value)/max_value*100
    # print(diff)
    ax2.plot(z.iloc[1,-1] , diff , '*',color=colors[k], lw=lwh)

ax2.set_ylabel('diff(%)',fontsize=12) 
ax2.set_ylim([0, 20])
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax2.legend(loc=4 , prop={'size': 10}) # 
# ax2.yaxis.label.set_color('red')
axes.set_xlabel('$\Gamma_1$',fontsize=12)
axes.set_ylabel('$P_2/P_1$',fontsize=12) 
axes.set_ylim([0, 2.5])
axes.set_title('$P_2/P_1$ under conditions sharing same $\Gamma_1$ with different $P_1,T_1$')
fig1.savefig("files/mm_similar_gp.eps")

# fig2 = plt.figure( dpi=300)
# lwh = 2
# axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# ax2 = axes.twinx()
# for k in range(0,len(Z1),1):
#     data = 'data' + str(k) + '.csv'
#     z = pd.read_csv(data, ",", skiprows=0)
#     axes.plot(z.iloc[:,-1] , z.iloc[:,4] , 'o',color=colors[k], lw=lwh)
#     max_value = z.iloc[:,4][np.argmax(abs(z.iloc[:,4]))]
#     min_value = z.iloc[:,4][np.argmin(abs(z.iloc[:,4]))]
#     diff = (max_value-min_value)/max_value*100
#     # print(diff)
#     ax2.plot(z.iloc[1,-1] , diff , '*',color=colors[k], lw=lwh)

# ax2.set_ylabel('diff(%)',fontsize=12) 
# # ax2.set_ylim([0, 10])
# ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# # ax2.legend(loc=4 , prop={'size': 10}) # 
# # ax2.yaxis.label.set_color('red')
# axes.set_xlabel('$Z_1$',fontsize=12)
# axes.set_ylabel('$T_2/T_1$',fontsize=12) 
# # axes.set_ylim([0.9, 1.1])
# axes.set_title('$T_2/T_1$ under conditions sharing same $Z_1$ with different $P_1,T_1$')
# # fig2.savefig("files/mm_similar_zt.eps")

fig3 = plt.figure( dpi=300)
lwh = 2
axes = fig3.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax2 = axes.twinx()
for k in range(0,len(G1),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)
    axes.plot(z.iloc[:,-1] , z.iloc[:,6] , 'o',color=colors[k], lw=lwh)
    max_value = z.iloc[:,6][np.argmax(abs(z.iloc[:,6]))]
    min_value = z.iloc[:,6][np.argmin(abs(z.iloc[:,6]))]
    diff = (max_value-min_value)/max_value*100
    # print(diff)
    ax2.plot(z.iloc[1,-1] , diff , '*',color=colors[k], lw=lwh)

ax2.set_ylabel('diff(%)',fontsize=12) 
ax2.set_ylim([0, 30])
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax2.legend(loc=4 , prop={'size': 10}) # 
# ax2.yaxis.label.set_color('red')
axes.set_xlabel('$\Gamma_1$',fontsize=12)
axes.set_ylabel('$M_2/M_1$',fontsize=12) 
axes.set_ylim([0.1, 0.5])
axes.set_title('$M_2/M_1$ under conditions sharing same $\Gamma_1$ with different $P_1,T_1$')
fig3.savefig("files/mm_similar_gm.eps")

fig4 = plt.figure( dpi=300)
lwh = 2
axes = fig4.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
for k in range(0,len(G1),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)
    axes.plot(z.iloc[:,-1] , z.iloc[:,-3] , 'o',color=colors[k], lw=lwh)
axes.set_xlabel('$\Gamma_1$',fontsize=12)
axes.set_ylabel('diff(%)',fontsize=12) 


