
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cv2
import subprocess
import glob
import sys
import platform


plt.rcParams["font.family"] = "Arial"
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams['ytick.major.width'] = 2.0
plt.rcParams['xtick.major.width'] = 2.0
plt.rcParams['xtick.minor.width'] = 1.4
plt.rcParams['ytick.minor.width'] = 1.4
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["ytick.major.size"] = 8
plt.rcParams["xtick.minor.size"] = 6
plt.rcParams["ytick.minor.size"] = 6
plt.rcParams['font.size'] = 24
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['figure.figsize'] = 10, 6

cm = plt.cm.get_cmap('tab20')
ctab = list(mcolors.TABLEAU_COLORS.values())
ctab = ['b', 'firebrick', 'y', 'g']


read_file_names = 'RESULTS_cen_T_H_ev5_csv/RESULTS_sph_flame_TH18*.csv'
read_file_list = glob.glob(read_file_names)
print(read_file_list)

fig, ax1 = plt.subplots(nrows=1, ncols=1, sharex='col', figsize=(10,7))
ax2 = ax1.twinx()

for i, read_file in enumerate(read_file_list):

    label = read_file.split('TH')[1].split('.csv')[0]
    df_cen = pd.read_csv(read_file)
    print(df_cen)

    df_cen.plot(kind='line', ax=ax1, x='t_s', y='T_c',  lw=4, marker=None, c=ctab[i], legend=False)
    df_cen.plot(kind='line', ax=ax2, x='t_s', y='YH_c', lw=4, marker=None, c=ctab[i], legend=False)

ax1.set_xlabel('t [s]')
ax1.set_ylabel('$T_c$ [K]')
ax2.set_ylabel('$Y_{H,c}$')

ax1.set_xlim(0, 0.001)
ax1.set_ylim(100, 2700)
ax2.set_ylim(-1.5e-4, 1.5e-3)

fig.tight_layout()
fig.savefig('time_cen_T_YH.png')

plt.show()

fig.clf()
plt.close('all')


