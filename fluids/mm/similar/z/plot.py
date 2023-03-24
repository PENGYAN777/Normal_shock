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
Z1 = np.arange(0.5, 0.92, 0.05).tolist()
for k in range(0,len(Z1),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)

    

"""
plot 
"""
n = 10
colors = plt.cm.tab20(np.linspace(0, 1, n))

fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax2 = axes.twinx()
for k in range(0,len(Z1),1):
    data = 'data' + str(k) + '.csv'
    z = pd.read_csv(data, ",", skiprows=0)
    axes.plot(z.iloc[:,-1] , z.iloc[:,3] , 'o',color=colors[k], lw=lwh)
    max_value = z.iloc[:,3][np.argmax(abs(z.iloc[:,3]))]
    min_value = z.iloc[:,3][np.argmin(abs(z.iloc[:,3]))]
    diff = (max_value-min_value)/max_value*100
    print(diff)
    ax2.plot(z.iloc[1,-1] , diff , '*',color=colors[k], lw=lwh)

ax2.set_ylabel('diff(%)',fontsize=12) 
ax2.set_ylim([0, 20])
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax2.legend(loc=4 , prop={'size': 10}) # 
# ax2.yaxis.label.set_color('red')
axes.set_xlabel('$Z_1$',fontsize=12)
axes.set_ylabel('$P_2/P_1$',fontsize=12) 
axes.set_ylim([0, 3])
axes.set_title('$P_2/P_1$ under conditions sharing same $Z_1$ with different $P_1,T_1$')
# axes.legend(loc=0 , prop={'size': 10}) # 

# ax2 = axes.twinx()  
# ax2.plot(data.iloc[:,2] , data.iloc[:,4] , 'ro', lw=lwh, label="$T_2/T_1$")
# ax2.set_ylabel('$T_2/T_1$',fontsize=12) 
# ax2.legend(loc=5 , prop={'size': 10}) # 
# ax2.yaxis.label.set_color('red')
# axes.set_title('$P_2/P_1$ vs $Z_t$',fontsize=14)
# axes.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
fig1.savefig("files/mm_similar_z.eps")



