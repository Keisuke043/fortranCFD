
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import glob
import os
import math
import subprocess

import settings_matplotlib as setMPL
import Savitzky_Golay_filter as SG


def calc_flamespeed_sg():

    time = np.array(df_peaks.loc[:,'time_conv'])
    xloc = np.array(df_peaks.loc[:,'xloc_conv'])

    window = 5
    xloc_sg = SG.non_uniform_savgol_der(time, xloc, window, 2, 0)

    if every_num == 10:
        window = 11
    if every_num == 5:
        window = 21
    if every_num == 1:
        window = 201
        window = 251
        # window = 321
    if len(time) < window and len(time) % 2 == 1:
        window = len(time)
    elif len(time) < window and len(time) % 2 == 0:
        window = len(time)-1
    
    dxdt_sg = SG.non_uniform_savgol_der(time, xloc, window, 2, 1)

    df_peaks.loc[:,'sb_cms_sg'] = dxdt_sg


def valid_convolve(xx, size):

    b = np.ones(size)/size
    xx_mean = np.convolve(xx, b, mode="same")

    n_conv = math.ceil(size/2)

    xx_mean[0] *= size/n_conv
    for i in range(1, n_conv):
        xx_mean[i] *= size/(i+n_conv)
        xx_mean[-i] *= size/(i + n_conv - (size % 2))

    return xx_mean


def calc_convolve():

    print(df_peaks)
    time = np.array(df_peaks.loc[:,'t_s'])
    xloc = np.array(df_peaks.loc[:,'xcm_1'])

    if every_num == 10:
        num = 10
    if every_num == 5:
        num = 20
    if every_num == 1:
        num = 30
        num = 50
    num = 1
        # num = 100

    b=np.ones(num)/num

    time_comv = valid_convolve(time, num)
    xloc_comv = valid_convolve(xloc, num)

    df_peaks['time_conv'] = time_comv
    df_peaks['xloc_conv'] = xloc_comv


def plot_r_t(save_dir):

    lw = 3
    mk = None
    if flag_marker:
        lw = 1
        mk = '.'

    fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(9,7))

    for i, (ener, df_peaks) in enumerate(df_peaks_dic.items()):
        e_label = str(ener)+' mJ'
        df_peaks.plot(kind='line', ax=ax, x='tms', y='xloc_conv', lw=lw,
                      marker=mk, ls='-', c=colors[i], legend=False, label=e_label)

    ax.set_xlim(tlim)
    ax.set_ylim(rlim)

    ax.set_xlabel('$t$ [ms]')
    ax.set_ylabel('$r$ [cm]')
    ax.legend(fontsize=14, loc='lower right')

    save_path = '{0}/r_t.png'.format(save_dir)

    fig.tight_layout()
    fig.savefig(save_path)
    fig.clf()

    subprocess.call(['imgcat',save_path])

    plt.close('all')


def plot_r_sb(save_dir):

    lw = 3
    mk = None
    if flag_marker:
        lw = 1
        mk = '.'

    fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(9,7))

    for i, (ener, df_peaks) in enumerate(df_peaks_dic.items()):
        e_label = str(ener)+' mJ'
        df_peaks.plot(kind='line', ax=ax, x='xloc_conv', y='sb_cms_sg', lw=lw,
                      marker=mk, ls='-', c=colors[i], legend=False, label=e_label)

    ax.set_xlim(rlim)
    ax.set_ylim(sblim)

    ax.set_xlabel('$r$ [cm]')
    ax.set_ylabel('$S_b$ [cm/s]')
    ax.legend(fontsize=14, loc='lower right')

    save_path = '{0}/r_sb.png'.format(save_dir)

    fig.tight_layout()
    fig.savefig(save_path)
    fig.clf()

    subprocess.call(['imgcat',save_path])

    plt.close('all')


def plot_k_sb(save_dir):

    lw = 3
    mk = None
    if flag_marker:
        lw = 1
        mk = '.'

    fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(9,7))

    for i, (ener, df_peaks) in enumerate(df_peaks_dic.items()):
        e_label = str(ener)+' mJ'
        df_peaks.plot(kind='line', ax=ax, x='kappa', y='sb_cms_sg', lw=lw,
                      marker=mk, ls='-', c=colors[i], legend=False, label=e_label)

    ax.set_xlim(klim)
    ax.set_ylim(sblim)

    ax.set_xlabel('$\kappa$ [1/s]')
    ax.set_ylabel('$S_b$ [cm/s]')
    ax.legend(fontsize=14, loc='lower right')

    save_path = '{0}/k_sb.png'.format(save_dir)

    fig.tight_layout()
    fig.savefig(save_path)
    fig.clf()

    subprocess.call(['imgcat',save_path])

    plt.close('all')


def plot_r_sb_alpha(save_dir):

    lw = 3
    mk = None
    if flag_marker:
        lw = 1
        mk = '.'
    alpha = 0.5

    for ener_solid in df_peaks_dic.keys():

        fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(9,7))

        for i, (ener, df_peaks) in enumerate(df_peaks_dic.items()):

            if ener == ener_solid:
                alpha = 1
            else:
                alpha = 0.2
            e_label = str(ener)+' mJ'
            df_peaks.plot(kind='line', ax=ax, x='xloc_conv', y='sb_cms_sg', lw=lw, alpha=alpha,
                          marker=mk, ls='-', c=colors[i], legend=False, label=e_label)

        ax.set_xlim(rlim)
        ax.set_ylim(sblim)

        ax.set_xlabel('$r$ [cm]')
        ax.set_ylabel('$S_b$ [cm/s]')
        ax.legend(fontsize=14, loc='lower right')

        save_path = '{0}/r_sb_{1}mJ.png'.format(save_dir, ener_solid)

        fig.tight_layout()
        fig.savefig(save_path)
        fig.clf()

        # subprocess.call(['imgcat',save_path])

    plt.close('all')


tab_colors_dark  = setMPL.tab_colors['dark']
tab_colors_light = setMPL.tab_colors['light']
colors = list(tab_colors_dark.values()) + list(tab_colors_light.values())
print(colors)


if __name__ == "__main__":

    flag_log_HRR = False
    flag_log_HRR = True
    flag_marker = True
    flag_marker = False
    # xlim = [3.1, 6.9] # cm
    # tlim = [0, None] # sec

    klim = [0.0, 1e+4] # 1/s
    rlim = [0.0, 0.8] # cm
    tlim = [0.0, 8.0] # msec
    sblim = [0.0, 350.0] # cm/s
    threshold_HRR = 1e+0 # J/m^3-s
    # klim = None # 1/s
    # rlim = None # cm
    # tlim = None # msec

    peak_param = 'HRR'
    extract_params = 'all'
    extract_params = ''

    every_num = 1
    every_num = 5

    tig_drop = 0.5  # ms  if sb < 0:
    tig = 15  # ms
    dir_base1 = 'RESULTS_sph_flame_peaks_ev{}_csv'.format(every_num)
    dir_base2 = 'isoOctane_SIP_phi0.8_p1_T300_rig500.0_tig200.0_*mJ_dx2.0d-5'
    dir_base2 = 'RESULTS_sph_flame_TH*'

    if extract_params == 'all':
        dir_base1 = dir_base1 + '_all'

    save_dir_base = 'RESULTS_sph_flame_peaks_ev{}_png'.format(every_num)
    subprocess.call(['mkdir','-p', save_dir_base])

    peak_num_list = [1, 2, 3]
    peak_num_list = [1]

    read_dirs = glob.glob(dir_base1+'/'+dir_base2)

    df_peaks_dic = {}

    read_dirs = sorted(read_dirs, reverse=True)
    print(read_dirs)

    for read_dir in read_dirs:

        dir_name = read_dir.split('/')[-1]
        temp_K  = float(dir_name.split('RESULTS_sph_flame_TH')[1])
        print(temp_K)

        # temp_K  = float(dir_name.split('_T')[1].split('_rig')[0])
        # ener_mJ = float(dir_name.split('_tig200.0_')[1].split('mJ')[0])
        # if ener_mJ == 0.2 or ener_mJ == 0.8:
        #     continue
        # if ener_mJ == 400.0 or ener_mJ == 800.0:
        #     continue

        save_dir = './{0}/{1}'.format(save_dir_base, dir_name)
        subprocess.call(['mkdir','-p', save_dir])

        print('read directory:', read_dir)
        print('save directory:', save_dir)

        read_path = '{0}/peaks_{1}.csv'.format(read_dir, peak_param)
        df_peaks = pd.read_csv(read_path)
        columns_bool = df_peaks.columns.str.contains('x_m_')
        columns_xm = df_peaks.loc[:, columns_bool].columns
        for col_xm in columns_xm:
            peak_num = col_xm.split('x_m_')[1]
            columns_name = 'xcm_{}'.format(peak_num)
            df_peaks[columns_name] = 1e+2*df_peaks[col_xm]
            df_peaks['tms'] = 1e+3*df_peaks['t_s']

        dt = df_peaks['t_s'].diff()[1]

        df_peaks['sb_cms'] = df_peaks['xcm_1'].diff()/dt

        idx_drop = df_peaks[df_peaks['sb_cms'] < 0].index.max()
        idx_drop = df_peaks[(df_peaks['tms'] < tig_drop) & (df_peaks['sb_cms'] < 0)].index.max()
        print('temp_K', temp_K, ', idx_drop', idx_drop)
        if not np.isnan(idx_drop):
            df_peaks = df_peaks[idx_drop < df_peaks.index]
        else:
            df_peaks = df_peaks.dropna(subset=['xcm_1'])

        for index in df_peaks.index:
            if df_peaks.loc[index, 'sb_cms'] == 0.0:
                df_peaks.drop(index, axis=0, inplace=True)
            else:
                break

        df_peaks = df_peaks[df_peaks['t_s']<tig*1e-3]

        calc_convolve()

        # df_peaks['sb_cms'] = df_peaks['xloc_conv'].diff()/dt
        # idx_drop = df_peaks[df_peaks['sb_cms'] <= 0].index.max()
        # if not np.isnan(idx_drop):
        #     df_peaks = df_peaks[idx_drop < df_peaks.index]

        # if ener_mJ <= 2.0:
        #     idx_drop = df_peaks[df_peaks['xcm_1'] <= 0.001].index.min()
        #     if not np.isnan(idx_drop):
        #         df_peaks = df_peaks[df_peaks.index < idx_drop]

        if df_peaks.empty:
            df_peaks['sb_cms_sg'] = np.nan
        else:
            calc_flamespeed_sg()

        df_peaks['kappa'] = 2.0*df_peaks['sb_cms_sg']/df_peaks['xloc_conv']
        # df_peaks_dic[ener_mJ] = df_peaks
        df_peaks_dic[temp_K] = df_peaks

        df_write = df_peaks.rename(columns = {'tms':'t [s]', 'xloc_conv':'r [cm]',
                                              'sb_cms_sg':'sb [cm/s]', 'kappa':'K [1/s]'})
        df_write = df_write[['t [s]','r [cm]','sb [cm/s]','K [1/s]']]
        df_write.to_csv('{}.csv'.format(save_dir), index=False)


    df_peaks_dic = dict(sorted(df_peaks_dic.items(), reverse=True))
    plot_r_t(save_dir_base)
    # plot_r_sb(save_dir_base)
    # plot_k_sb(save_dir_base)
    # plot_r_sb_alpha(save_dir_base)



