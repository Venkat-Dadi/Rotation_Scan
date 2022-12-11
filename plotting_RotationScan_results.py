#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 17:37:01 2022

@author: venkatareddy
"""


import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.image as mpimage
import matplotlib.lines as lines
import matplotlib.colors as colors
import sys
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterSciNotation
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, FixedLocator,LogFormatter )
matplotlib.rcParams['axes.linewidth'] = 0.5



# =============================================================================
# if len(sys.argv) == 1:
#     print('Ploting Output from ProteinEuler program.')
#     print('Usage: python plot_results_ProteinEuler_4.py results.data') 
#     print('\n')
#     sys.exit()
# =============================================================================

input_data = sys.argv[1]
#input_data = 'results_PE_dot1L-6nj9-randOri1_NCP-6nj9_res010_20220209_15h30m44s.data'
f1 = open(input_data, 'r')

data1 = []
for line in f1.readlines():
    if not line.startswith('#'):
        data1.append(line.replace('\t', ',').replace(' ', '').strip('\n'))

data2 = []        
for d1 in data1:
    data2.append(np.fromstring(d1, sep=','))
    
data3 = np.array(data2)

header = ['alpha', 'beta', 'gamma', 'ENZ_CA', 'SUB_CA', 'TOT_CA', 'clashes'] 

df1 = pd.DataFrame(data3.astype(dtype=int), columns=header)

df2 = df1.sort_values(by=['TOT_CA'], ascending=True)

output = input_data.split('/')[-1].strip('\.data')
df2.to_csv('{}_sorted.csv'.format(output), index=False)

min_alpha = df2['alpha'][:1].tolist()[0]
min_beta = df2['beta'][:1].tolist()[0]
min_gamma = df2['gamma'][:1].tolist()[0]


min_alpha_df = df2[df2['alpha'] == int(min_alpha)]
min_beta_df = df2[df2['beta'] == int(min_beta)]
min_gamma_df = df2[df2['gamma'] == int(min_gamma)]

min_alpha_df_TOTCA = min_alpha_df['TOT_CA'].tolist()
min_beta_df_TOTCA = min_beta_df['TOT_CA'].tolist()
min_gamma_df_TOTCA = min_gamma_df['TOT_CA'].tolist()

min_TOTCA = min([min(min_alpha_df_TOTCA), min(min_beta_df_TOTCA), min(min_gamma_df_TOTCA)])
max_TOTCA = max([max(min_alpha_df_TOTCA), max(min_beta_df_TOTCA), max(min_gamma_df_TOTCA)])

min_alpha_df2 = pd.DataFrame(min_alpha_df, columns=['beta', 'gamma', 'TOT_CA'])
min_alpha_df2_pivot = min_alpha_df2.pivot('beta', 'gamma', 'TOT_CA')

min_beta_df2 = pd.DataFrame(min_beta_df, columns=['alpha', 'gamma', 'TOT_CA'])
min_beta_df2_pivot = min_beta_df2.pivot('alpha', 'gamma', 'TOT_CA')

min_gamma_df2 = pd.DataFrame(min_gamma_df, columns=['alpha', 'beta', 'TOT_CA'])
min_gamma_df2_pivot = min_gamma_df2.pivot('alpha', 'beta', 'TOT_CA')


l = list(set(df2['alpha'].tolist()))
l.sort()
res = abs(l[0]-l[1])


alpha_range = np.arange(df2['alpha'].min(), df2['alpha'].max()+1, res).tolist()
beta_range = np.arange(df2['beta'].min(), df2['beta'].max()+1, res).tolist()
gamma_range = np.arange(df2['gamma'].min(), df2['gamma'].max()+1, res).tolist()

alpha = 'alpha'
beta = 'beta'
gamma = 'gamma'

#fig = plt.figure(figsize=(7.086617, 4), dpi = 600)
fig = plt.figure(figsize=(7.4, 4.4))

ax1 = fig.add_axes([0.08, 0.17, 0.2, 0.7])
ax2 = fig.add_axes([0.30, 0.17, 0.4, 0.7])
ax3 = fig.add_axes([0.72, 0.17, 0.2, 0.7])
cbaxes1 = fig.add_axes([0.08, 0.12, 0.84, 0.03])  

font_size = 8

orig_cmap = matplotlib.cm.bwr

ax1.set_xticks(np.arange(len(beta_range)))
ax1.set_yticks(np.arange(len(alpha_range)))
ax1.set_xticklabels(beta_range, fontsize = font_size, rotation=90)
ax1.set_yticklabels(alpha_range, fontsize = font_size)
ax1.set_xlabel(r'$\beta~(^\circ)$', fontsize = font_size)
ax1.set_ylabel(r'$\gamma~(^\circ)$', fontsize = font_size)
im_1 = ax1.imshow(min_alpha_df2_pivot.T, cmap=orig_cmap, aspect = 'auto',\
                  norm=colors.LogNorm())
ax1.invert_yaxis()
#ax1.set_title('(A)')
ax1.text(6, 34, r'$\alpha= {}^\circ$'.format(min_alpha), fontsize = font_size)

ax2.set_xticks(np.arange(len(alpha_range)))
ax2.set_yticks(np.arange(len(gamma_range)))
ax2.set_yticklabels([])
ax2.set_xticklabels(alpha_range, fontsize = font_size, rotation=90)
#ax2.yaxis.tick_right()
#ax2.set_yticklabels(gamma_range, fontsize = 12)
ax2.set_xlabel(r'$\alpha~(^\circ)$', fontsize = font_size)
#ax2.set_ylabel(r'$\beta~(^\circ)$', fontsize = 14)
im_2 = ax2.imshow(min_beta_df2_pivot.T, cmap=orig_cmap, aspect = 'auto',\
                  norm=colors.LogNorm())
ax2.invert_yaxis()
ax2.text(15, 34, r'$\beta= {}^\circ$'.format(min_beta), fontsize = font_size)

ax3.set_xticks(np.arange(len(beta_range)))
ax3.set_yticks(np.arange(len(gamma_range)))
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position('right')
#ax3.set_yticks([])
ax3.set_xticklabels(beta_range, fontsize = font_size, rotation=90)
ax3.set_yticklabels(gamma_range, fontsize = font_size)
ax3.set_xlabel(r'$\beta~(^\circ)$', fontsize = font_size)
ax3.set_ylabel(r'$\alpha~(^\circ)$', fontsize = font_size)
im_3 = ax3.imshow(min_gamma_df2_pivot, cmap=orig_cmap, aspect = 'auto',\
                  norm=colors.LogNorm())
ax3.invert_yaxis()
ax3.text(6, 34, r'$\gamma= {}^\circ$'.format(min_gamma), fontsize = font_size)


ax = [ax1, ax2, ax3]

ticks_major = [ee for ee in np.arange(len(gamma_range)) if ee%2 == 0]
ticks_minor = [ee for ee in np.arange(len(gamma_range)) if ee%2 != 0]

for a3 in ax:
    a3.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(ticks_major))
    a3.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(ticks_minor))
    a3.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(ticks_major))
    a3.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(ticks_minor))

for c1 in ax:
    c1.xaxis.tick_top()
    
for l1 in ax:
    l1.xaxis.set_label_position('top')    
    
every_nth = 2
for a1 in ax:
    for n1, label1 in enumerate(a1.xaxis.get_ticklabels()):
        if n1 % every_nth != 0:
            label1.set_visible(True)
for a2 in ax:
    for n2, label2 in enumerate(a2.yaxis.get_ticklabels()):
        if n2 % every_nth != 0:
            label2.set_visible(True)   
            
for t1 in ax:
    t1.xaxis.set_tick_params(width=0.5)  
    t1.yaxis.set_tick_params(width=0.5)     

formatter = LogFormatter(10, labelOnlyBase=False) 
cbar_ticks = [min_TOTCA, 100, 1000, max_TOTCA]  
cbar = plt.colorbar(im_1, cax = cbaxes1, orientation='horizontal',\
                    ticks = cbar_ticks, format=formatter)
cbar.ax.set_xticklabels(cbar_ticks)
cbar.ax.tick_params(labelsize = font_size, width = 0.5)
#cbar.set_label(label = 'TCA', size = font_size, loc ='right')
#cbar.ax.set_title('TCA', fontsize=font_size)
cbar.ax.text(max_TOTCA+250, 0.20, 'TCA', fontsize=font_size)
#plt.show()
plt.savefig('{}.pdf'.format(output), facecolor='none', dpi=600, transparent=True)  
