
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import glob
import os
import subprocess

import math

import settings_matplotlib as setMPL
import Savitzky_Golay_filter as SG


def calc_flamespeed_sg():

    # time = np.array(df_peaks.loc[:,'t_s'])
    # xloc = np.array(df_peaks.loc[:,'xcm_1'])
    time = np.array(df_peaks.loc[:,'time_conv'])
    xloc = np.array(df_peaks.loc[:,'xloc_conv'])

    window = 5
    xloc_sg = SG.non_uniform_savgol_der(time, xloc, window, 2, 0)
    window = 21
    window = 11
    dxdt_sg = SG.non_uniform_savgol_der(time, xloc_sg, window, 2, 1)

    df_peaks.loc[:,'xcm_sg'] = xloc_sg
    df_peaks.loc[:,'sb_cms_sg'] = dxdt_sg


def valid_convolve(xx, size):

    b = np.ones(size)/size
    xx_mean = np.convolve(xx, b, mode="same")

    n_conv = math.ceil(size/2)

    # 補正部分
    xx_mean[0] *= size/n_conv
    for i in range(1, n_conv):
        xx_mean[i] *= size/(i+n_conv)
        xx_mean[-i] *= size/(i + n_conv - (size % 2))
	# size%2は奇数偶数での違いに対応するため

    return xx_mean


def calc_convolve():

    time = np.array(df_peaks.loc[:,'t_s'])
    xloc = np.array(df_peaks.loc[:,'xcm_1'])
    # sbcms = np.array(df_peaks['sb_cms_sg'])
    # x = np.array(df_peaks.loc[:,'t_s'])
    # y = np.array(df_peaks.loc[:,'xcm_1'])

    # x = df_peaks['xcm_sg']
    # y = df_peaks['sb_cms_sg']
    # print(len(x), len(y))

    num = 100
    num = 200
    num = 800
    num = 30
    num = 5
    num = 10

    b=np.ones(num)/num

    # xloc_comv=np.convolve(xloc, b, mode='same')

    time_comv = valid_convolve(time, num)
    xloc_comv = valid_convolve(xloc, num)

    df_peaks['time_conv'] = time_comv
    df_peaks['xloc_conv'] = xloc_comv
    # df_peaks['xloc_conv'] = xloc_comv
    # df_peaks['sb_cms_sg_conv'] = sbcms_comv
    # df_peaks['xloc'] = xloc_comv
    # df_peaks['sb_cms_sg'] = sbcms_comv

    # print(len(y2))
    # #移動平
    # #y2 = inputData.rolling('x_m_1',20,center=True).sum()
    # #F = np.fft.fft(y)

    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # ax1.yaxis.set_ticks_position('both')
    # ax1.xaxis.set_ticks_position('both')
    # ax1.set_xlabel('Time [ms]')
    # ax1.set_ylabel('radius')
    # ax1.plot(x, y, label="original")
    # ax1.plot(x, y2, label="convolve")
    # 
    # 
    # fig.tight_layout()
    # plt.legend(loc='upper left')
    # # plt.savefig(dataName+'_convolve.png')

    # plt.show()
    # plt.close()
    # 
    # # df = pd.DataFrame({'x': x, 'y': y2})
    # # df.to_csv(dataName+'_convolve.csv')


def graph_t_x(save_dir):

    fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(9,7))

    lw = 3
    mk = None
    if flag_marker:
        lw = 1
        mk = '.'

    colors = ['r', 'b', 'g', 'y']
    labels = ['1st', '2nd', '3rd', '4th']

    for peak_num in reversed(peak_num_list):
        x_columns = 'xcm_{}'.format(peak_num)
        label = labels[peak_num-1]
        color = colors[peak_num-1]
        df_peaks.plot(kind='line', ax=ax, x='tms', y=x_columns, lw=lw,
                 marker=mk, ls='-', c=color, legend=False, label=label)

    ax.set_xlim(tlim)
    ax.set_ylim(rlim)

    ax.set_xlabel('t [ms]')
    ax.set_ylabel('r [cm]')

    save_path = '{0}/t_x_peaknum{1}.png'.format(save_dir, peak_num_list[-1])

    fig.tight_layout()
    fig.savefig(save_path)
    fig.clf()

    plt.close('all')


def plot_x_t(save_dir):

    lw = 3
    mk = None
    if flag_marker:
        lw = 1
        mk = '.'

    fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(9,7))

    for i, (ener, df_peaks) in enumerate(df_peaks_dic.items()):
        e_label = str(ener)+' mJ'
        df_peaks.plot(kind='line', ax=ax, x='tms', y='xcm_sg', lw=lw,
                      marker=mk, ls='-', c=colors[i], legend=False, label=e_label)

    ax.set_xlim(tlim)
    ax.set_ylim(rlim)

    ax.set_xlabel('t [ms]')
    ax.set_ylabel('r [cm]')
    ax.legend(fontsize=16, loc='upper right')

    save_path = '{0}/x_t.png'.format(save_dir)

    fig.tight_layout()
    fig.savefig(save_path)
    fig.clf()

    subprocess.call(['imgcat',save_path])

    plt.close('all')


def plot_x_sb(save_dir):

    lw = 3
    mk = None
    if flag_marker:
        lw = 1
        mk = '.'

    fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(9,7))

    for i, (ener, df_peaks) in enumerate(df_peaks_dic.items()):
        e_label = str(ener)+' mJ'
        df_peaks.plot(kind='line', ax=ax, x='xcm_sg', y='sb_cms_sg', lw=lw,
                      marker=mk, ls='-', c=colors[i], legend=False, label=e_label)

    ax.set_xlim(rlim)
    ax.set_ylim(sblim)

    ax.set_xlabel('r [cm]')
    ax.set_ylabel('Sb [cm/s]')
    ax.legend(fontsize=18, loc='lower right')

    save_path = '{0}/x_sb.png'.format(save_dir)

    fig.tight_layout()
    fig.savefig(save_path)
    fig.clf()

    subprocess.call(['imgcat',save_path])

    plt.close('all')


tab_colors_dark  = setMPL.tab_colors['dark']
tab_colors_light = setMPL.tab_colors['light']
colors = list(tab_colors_dark.values()) + list(tab_colors_light.values())
print(colors)


if __name__ == "__main__":

    flag_log_HRR = False
    flag_log_HRR = True
    flag_marker = False
    flag_marker = True
    # xlim = [3.1, 6.9] # cm
    # tlim = [0, None] # sec
    # rlim = [0.0, 4.0] # cm
    rlim = [0.0, 0.5] # cm
    rlim = [0.0, 0.4] # cm
    tlim = [0.0, 3.0] # msec
    sblim = [0.0, 500.0] # cm/s
    threshold_HRR = 1e+0 # J/m^3-s

    peak_param = 'HRR'
    extract_params = 'all'
    extract_params = ''

    dir_base1 = 'RESULTS_sph_flame_peaks_csv'
    dir_base2 = 'RESULTS_sph_flame_*mJ'

    if extract_params == 'all':
        dir_base1 = dir_base1 + '_all'

    save_dir_base = 'RESULTS_sph_flame_peaks_png'.format(peak_param)

    # graph_func = graph_contour
    graph_func = graph_t_x

    peak_num_list = [1, 2]
    peak_num_list = [1, 2, 3]
    peak_num_list = [1]

    read_dirs = glob.glob(dir_base1+'/'+dir_base2)

    df_peaks_dic = {}

    read_dirs = sorted(read_dirs, reverse=True)

    for read_dir in read_dirs:

        dir_name = read_dir.split('/')[-1]
        ener_mJ = float(dir_name.split('sph_flame_')[1].split('mJ')[0])
        if ener_mJ < 0.145:
            continue
        # if not ener_mJ == 0.1:
        #     continue

        # print(ener_mJ)

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
        if not np.isnan(idx_drop):
            df_peaks = df_peaks[idx_drop < df_peaks.index]

        for index in df_peaks.index:
            if df_peaks.loc[index, 'sb_cms'] == 0.0:
                df_peaks.drop(index, axis=0, inplace=True)
            else:
                break

        calc_convolve()

        # df_peaks['sb_cms'] = df_peaks['xloc_conv'].diff()/dt
        # idx_drop = df_peaks[df_peaks['sb_cms'] <= 0].index.max()
        # print(idx_drop)
        # if not np.isnan(idx_drop):
        #     df_peaks = df_peaks[idx_drop < df_peaks.index]

        calc_flamespeed_sg()

        graph_func(save_dir)

        df_peaks_dic[ener_mJ] = df_peaks


        # df_dict = {}
        # df_all = pd.DataFrame()
        # for peak_num in reversed(peak_num_list):
        #     csv_file = './{0}/peaks_{1}_{2}.csv'.format(read_dir, peak_param, peak_num)
        #     df = pd.read_csv(csv_file)
        #     columns_m = 'x_m_{}'.format(peak_num)
        #     columns_cm = 'xcm_{}'.format(peak_num)
        #     df[columns_cm] = 1e+2*df[columns_m]
        #     df['tms'] = 1e+3*df['t_s']
        #     columns_list = ['x_m', 'xcm', 'HRR']
        #     for columns in columns_list:
        #         df = df.rename(columns={'{}_{}'.format(columns, peak_num):columns})
        #     df_all = pd.concat([df_all, df])
        # df_all = df_all[threshold_HRR < df_all['HRR']]

    plot_x_t(save_dir_base)
    plot_x_sb(save_dir_base)

