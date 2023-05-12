
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

import settings_matplotlib as setMPL


cm = plt.cm.get_cmap('tab20')
ctab = list(mcolors.TABLEAU_COLORS.values())


def pTrhou(save_dir, ev_num=1):

    for i, df in data_nout_dict.items():

        fig, ax = plt.subplots(nrows=3,ncols=1, sharex='col', figsize=(12,12))
        ax_0t = ax[0].twinx()
        ax_1t = ax[1].twinx()

        lw = 4
        mk = None
        if flag_marker:
            lw = 1
            mk = '.'

        df.plot(kind='line', ax=ax[0], x='x_m', y='  T', lw=lw, marker=mk, c=ctab[0], legend=False, label='T')
        df.plot(kind='line', ax=ax_0t, x='x_m', y='  p', lw=lw, marker=mk, c=ctab[1], legend=False, label='p')
        df.plot(kind='line', ax=ax[1], x='x_m', y='u_x', lw=lw, marker=mk, c=ctab[2], legend=False, label='u')
        df.plot(kind='line', ax=ax_1t, x='x_m', y='rho', lw=lw, marker=mk, c=ctab[3], legend=False, label='rho')
        ax[2].plot(df['x_m'], df['rho']*df['u_x'], lw=lw, marker=mk, c=ctab[4], label='rho*u')

        ax[0].set_ylabel('T [K]')
        ax_0t.set_ylabel('p [Pa]')
        ax[1].set_ylabel('u [m/s]')
        ax_1t.set_ylabel('rho [kg/m3]')
        ax[2].set_ylabel('rhou [kg/m$^2$-s]')
        ax[2].set_xlabel('x [m]')

        ax[0].set_ylim(Tlim)
        ax_0t.set_ylim(plim)
        ax[1].set_ylim(ulim)
        ax_1t.set_ylim(rholim)
        ax[2].set_ylim(rhoulim)
        ax[2].set_xlim(xlim)
       
        handler_0o, label_0o = ax[0].get_legend_handles_labels()
        handler_0t, label_0t = ax_0t.get_legend_handles_labels()
        ax[0].legend(handler_0o+handler_0t, label_0o+label_0t, loc='upper right', fontsize=16)
        handler_1o, label_1o = ax[1].get_legend_handles_labels()
        handler_1t, label_1t = ax_1t.get_legend_handles_labels()
        ax[1].legend(handler_1o+handler_1t, label_1o+label_1t, loc='upper right', fontsize=16)
        ax[2].legend(loc='upper right', fontsize=16)

        fig.suptitle('t = {:1.5f} s'.format(float(attr_nout_dict[i]['time from 0s'])), fontsize=18)
        fig.tight_layout()
        fig.savefig('{0}/nout_{1:05}.png'.format(save_dir, num_sta+i))
        fig.clf()
    plt.close('all')


def get_filenum():

    num_list = []
    filepath_num_dic = {}
    for i in range(0, len(file_list)):
        num_list.append(int(file_list[i].split('/')[-1].split('_nout')[1].split('.h5')[0]))
        filepath_num_dic[num_list[i]] = file_list[i]

    max_num_png_list = []
    if len(sys.argv)==2:
        subprocess.call(['rm','-r', save_dir])
    subprocess.call(['mkdir','-p', save_dir])
    file_png_list=glob.glob('{}/nout*'.format(save_dir))
    num_png_list = []
    if len(file_png_list) != 0:
        for i in range(0, len(file_png_list)):
            num_png_list.append(int(file_png_list[i].split('/')[-1].split('nout_')[1].split('.png')[0]))
        max_num_png_list.append(max(num_png_list))
    else:
        max_num_png_list.append(0)

    num_list = sorted(num_list)
    num_list_pop = [i for i in num_list if i <= min(max_num_png_list)*every_num_graph]
    idx_numlistpop_max = num_list.index(max(num_list_pop))
    num_list = num_list[idx_numlistpop_max:]
    print(num_list)

    return num_list, filepath_num_dic


def read_HDFfiles(i_data):

    data_nout_dict = {}
    attr_nout_dict = {}

    print('reading result files:', filepath_dic[i])
    f = h5py.File(filepath_dic[i], 'r')
    for nout_group in f.keys():
        df = pd.DataFrame()
        sr = pd.DataFrame()
        number_out =int(nout_group.split('_')[0].split('nout')[1])
        if (number_out%every_num_graph == 0):
            print('file_nout, fig_nout, numfig = {}, {}, {}'.format(i, number_out, num_sta+i_data))
            for key_attr in f[nout_group].attrs.keys():
                # print(key_attr)
                if key_attr == 'time from 0s':
                    sr[key_attr] = f[nout_group].attrs[key_attr]
                if key_attr == 'cpu time total':
                    sr[key_attr] = f[nout_group].attrs[key_attr]
                if key_attr == 'dt':
                    sr[key_attr] = f[nout_group].attrs[key_attr]
            attr_nout_dict[i_data]=sr
            for key_dataset in f[nout_group].keys():
                # print(key_dataset)
                group_path = f[nout_group+'/'+key_dataset].name
                sr_tmp = pd.Series(np.array(f[group_path]), name=key_dataset)
                df_tmp = pd.DataFrame(sr_tmp)
                df = pd.concat([df, df_tmp], axis=1)
            data_nout_dict[i_data]=df
            i_data += 1

    return i_data, data_nout_dict, attr_nout_dict


def mkmovie(save_dir):

    file_mp4 = '{0}.mp4'.format(save_dir)
    pf = platform.system()

    if pf == 'Darwin':

        cmd='ffmpeg -y -start_number 0 -r 30 -i {0}/nout_%05d.png -vcodec libx264 -pix_fmt yuv420p {1}'\
            .format(save_dir, file_mp4)
        subprocess.call(cmd.split())

    elif pf == 'Linux':

        frame_rate=30
        filelist = glob.glob('{0}/nout_*.png'.format(save_dir))
        filelist.sort()
        img = cv2.imread(filelist[0])
        height, width, color = img.shape
        fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
        video = cv2.VideoWriter(file_mp4, fourcc, frame_rate, (width, height))
        for file in filelist:
            img = cv2.imread(file)
            video.write(img)


if __name__ == "__main__":

    # print(sys.argv)
    
    ### setting params
    read_dir_name = 'RESULTS_SOD_hdf'
    flag_marker = False
    every_num_graph = 5
    xlim = [0, 1.0]
    Tlim = [200, 500]
    plim = [4e+4, 1.2e+5]
    ulim = [-10, 100]
    rholim = [0, 1.3]
    rhoulim = [-10, 90]


    read_path = '{}'.format(read_dir_name)
    file_list=glob.glob('{}/results*'.format(read_path))
    print(file_list)

    dir_name = read_dir_name.split('/')[-1].split('_hdf')[0]
    save_dir_base = dir_name+'_png'
    save_dir = './{0}/{1}_pTrhou_png'.format(save_dir_base, dir_name)
    print('save directory:', save_dir_base)

    num_list, filepath_dic = get_filenum()
    num_sta = int(num_list[0]/every_num_graph)

    i_data = 0
    for i in num_list:
        if i == num_list[-1]:
            try:
                f = h5py.File(filepath_dic[i], 'r')
            except:
                continue
        i_data, data_nout_dict, attr_nout_dict = read_HDFfiles(i_data)
        pTrhou(save_dir, ev_num=every_num_graph)
    mkmovie(save_dir)

