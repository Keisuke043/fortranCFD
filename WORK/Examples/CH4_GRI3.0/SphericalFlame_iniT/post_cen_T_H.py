
import os
import sys
import glob
import subprocess
import numpy as np
import pandas as pd
import cantera as ct
import h5py


def get_filenum(file_list):

    num_list = []
    filepath_num_dic = {}
    for i in range(0, len(file_list)):
        num_list.append(int(file_list[i].split('/')[-1].split('_nout')[1].split('.h5')[0]))
        filepath_num_dic[num_list[i]] = file_list[i]
    num_list = sorted(num_list)

    return num_list, filepath_num_dic


def read_HDFfiles(i_numfile, i_num_post, num_sta):

    data_nout_dict = {}
    attr_nout_dict = {}

    print('reading result files:', filepath_dic[i_numfile])
    f = h5py.File(filepath_dic[i_numfile], 'r')
    for nout_group in f.keys():
        df = pd.DataFrame()
        sr = pd.DataFrame()
        i_num =int(nout_group.split('_')[0].split('nout')[1])
        if (i_num%every_num == 0):
            print('file_nout, csv_nout, cur_num = {}, {}, {}'.format(i_numfile, i_num, num_sta+i_num_post))
            for key_attr in f[nout_group].attrs.keys():
                # print(key_attr)
                if key_attr == 'time from 0s':
                    sr[key_attr] = f[nout_group].attrs[key_attr]
                if key_attr == 'cpu time total':
                    sr[key_attr] = f[nout_group].attrs[key_attr]
                if key_attr == 'dt':
                    sr[key_attr] = f[nout_group].attrs[key_attr]
            attr_nout_dict[i_num_post]=sr
            for key_dataset in f[nout_group].keys():
                # print(key_dataset)
                group_path = f[nout_group+'/'+key_dataset].name
                sr_tmp = pd.Series(np.array(f[group_path]), name=key_dataset)
                df_tmp = pd.DataFrame(sr_tmp)
                df = pd.concat([df, df_tmp], axis=1)
            df = df.assign(xcm = 1e+2*df['x_m'])

            gas = ct.Solution(chem_cti)
            sp_names = gas.kinetics_species_names

            Yi_sp_name = ['Y_'+name for  name in gas.species_names]
            df_Yi = pd.DataFrame(columns=Yi_sp_name)

            hrr_ct_list = []
            for index, row in df.iterrows():
                T = row['  T']
                p = row['  p']
                Xi_dict = {sp:row[sp] for sp in sp_names}
                gas.TPX = T, p, Xi_dict

                hrr_ct = gas.heat_release_rate
                hrr_ct_list.append(hrr_ct)

                df_Yi.loc[index] = gas.Y

            df = df.assign(HRRct=hrr_ct_list)
            df = pd.concat([df, df_Yi], axis=1)

            data_nout_dict[i_num_post]=df
            i_num_post += 1

    return i_num_post, data_nout_dict, attr_nout_dict


def extract_T_H():

    # print(filenum_list)
    num_sta = 0
    i_data = 0
    df_cen = pd.DataFrame(index=[], columns=['t_s', 'T_c', 'XH_c', 'YH_c'])
    for i_numfile in filenum_list:
        if i_numfile == filenum_list[-1]:
            try:
                f = h5py.File(filepath_dic[i], 'r')
            except:
                continue
        i_data, data_nout_dict, attr_nout_dict = read_HDFfiles(i_numfile, i_data, num_sta)
        # print('Writing in csv file {}'.format(i_data))
        for i_nout, df in data_nout_dict.items():
            
            time = float(attr_nout_dict[i_nout]['time from 0s'])
            df_cen.loc[i_nout] = [time, df.loc[0, '  T'], df.loc[0, 'H'], df.loc[0, 'Y_H']]

    return df_cen


if __name__ == "__main__":

    read_path_name = 'RESULTS_sph_flame_TH*_hdf'
    read_path_list = glob.glob(read_path_name)
    print(read_path_list)

    chem_cti = './GRI3.0_CH4/chem.cti' # cti file for cantera postprocessing

    num_peaks = 3
    every_num = 10
    every_num = 5

    for read_path in read_path_list:

        file_list = glob.glob('{}/results_*.h5'.format(read_path))
        filenum_list, filepath_dic = get_filenum(file_list)
        f_hdf = h5py.File(filepath_dic[filenum_list[0]], 'r')

        print(read_path)
        print(filenum_list)

        save_dir = './RESULTS_cen_T_H_ev{}_csv'.format(every_num)
        subprocess.call(['mkdir','-p', save_dir])

        df_cen = extract_T_H()

        print(save_dir)
        print(df_cen)

        TH_dir = read_path.split('_hdf')[0]
        save_path = save_dir + '/{}.csv'.format(TH_dir)
        df_cen.to_csv(save_path, index=False)

        f_hdf.close()

