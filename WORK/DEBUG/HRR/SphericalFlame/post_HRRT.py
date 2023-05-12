
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
import cantera as ct


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


def graph_HRRT(save_dir, ev_num=1):

    for j, df in data_nout_dict.items():

        fig, ax1 = plt.subplots(nrows=1,ncols=1, figsize=(10,6))
        ax2 = ax1.twinx()

        mk = None
        if flag_marker:
            mk = '.'

        df.plot(ax=ax1, x='xcm', y='HRR', ls='--', lw=2.5, marker=mk, c='r', alpha=0.9, legend=False, label='HRR (Enthalpy diff)')
        df.plot(ax=ax1, x='xcm', y='HRD', ls='-.', lw=2.5, marker=mk, c='y', alpha=0.9, legend=False, label='HRR (Omega dot)')
        df.plot(ax=ax1, x='xcm', y='HRRct', ls='-', lw=2.5, marker=mk, c='c', alpha=0.9, legend=False, label='HRR (Cantera postprocessing)', zorder=-1)
        df.plot(ax=ax2, x='xcm', y='  T', ls='-', lw=2.5, marker=mk, c='k', alpha=0.9, legend=False)

        ax1.set_zorder(ax2.get_zorder()+1)
        ax1.patch.set_visible(False)

        ax1.set_xlabel('r [cm]')
        ax1.set_ylabel('HRR [J/m$^3$-s]')
        ax2.set_ylabel('T [K]')

        ax1.set_xlim(xlim)
        ax1.set_ylim(HRRlim)
        ax2.set_ylim(Tlim)

        if flag_log:
            ax1.set_ylim(HRRlim_log)
            ax1.set_yscale('log')

        ax1.set_title('t = {:1.5f} s'.format(float(attr_nout_dict[j]['time from 0s'])), fontsize=18)
        ax1.legend(loc='upper right', fontsize=13)

        fig.tight_layout()
        fig.savefig('{0}/nout_{1:05}.png'.format(save_dir, num_sta+j))
        fig.clf()
    plt.close('all')


def get_filenum():

    num_list = []
    filepath_num_dic = {}
    for i in range(0, len(file_list)):
        num_list.append(int(file_list[i].split('/')[-1].split('_nout')[1].split('.h5')[0]))
        filepath_num_dic[num_list[i]] = file_list[i]

    max_num_png_list = []
    if len(sys.argv) == 2:
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
        number_out = int(nout_group.split('_')[0].split('nout')[1])
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
            attr_nout_dict[i_data] = sr
            for key_dataset in f[nout_group].keys():
                # print(key_dataset)
                group_path = f[nout_group+'/'+key_dataset].name
                sr_tmp = pd.Series(np.array(f[group_path]), name=key_dataset)
                df_tmp = pd.DataFrame(sr_tmp)
                df = pd.concat([df, df_tmp], axis=1)
            df = df.assign(xcm=1e+2*df['x_m'])

            gas = ct.Solution(chem_cti)
            sp_names = gas.kinetics_species_names

            hrr_ct_list = []
            for index, row in df.iterrows():
                T = row['  T']
                p = row['  p']
                Xi_dict = {sp:row[sp] for sp in sp_names}
                gas.TPX = T, p, Xi_dict

                hrr_ct = gas.heat_release_rate
                hrr_ct_list.append(hrr_ct)
            df = df.assign(HRRct=hrr_ct_list)

            data_nout_dict[i_data]=df
            i_data += 1

    return i_data, data_nout_dict, attr_nout_dict


def mkmovie(save_dir):

    file_mp4 = '{0}.mp4'.format(save_dir)
    pf = platform.system()

    if pf == 'Darwin':

        cmd='ffmpeg -y -start_number 0 -r 30 -i {0}/nout_%05d.png -vcodec libx264 -pix_fmt yuv420p ./{1}'\
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
    read_dir_name = 'RESULTS_sph_flame_hdf'
    chem_cti = './Lusk39_DME/chem.cti' # cti file for postprocessing of cantera HRR
    flag_log = True
    flag_marker = False

    every_num_graph = 10

    xlim = [0, 0.5] # cm
    Tlim = [600, 3300] # K
    HRRlim = [0, 5e+10] # J/m^3-s
    HRRlim_log = [5e-1, 6e+12] # J/m^3-s

    read_path = '{}'.format(read_dir_name)
    file_list=glob.glob('{}/results*'.format(read_path))
    print(file_list)

    dir_name = read_dir_name.split('/')[-1]
    save_dir_base = dir_name+'_png'
    save_dir_sub = dir_name.split('RESULTS_')[1].split('_hdf')[0]
    save_dir_sub = save_dir_sub+'_HRRT'

    if flag_log:
        save_dir_sub = save_dir_sub+'_log'
    if flag_marker:
        save_dir_sub = save_dir_sub+'_marker'

    save_dir = './{0}/{1}_png'.format(save_dir_base, save_dir_sub)
    print('save directory:', save_dir)

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
        graph_HRRT(save_dir, ev_num=every_num_graph)
    mkmovie(save_dir)

