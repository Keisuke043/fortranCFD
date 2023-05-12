
import h5py
import pandas as pd
import numpy as np
import subprocess
import glob
import sys


def get_filenum():

    num_list = []
    filepath_num_dic = {}
    for i in range(0, len(file_list)):
        num_list.append(int(file_list[i].split('/')[-1].split('_nout')[1].split('.h5')[0]))
        filepath_num_dic[num_list[i]] = file_list[i]

    max_num_csv_list = []
    if len(sys.argv) == 2:
        subprocess.call(['rm','-r', save_dir])
    subprocess.call(['mkdir','-p', save_dir])
    file_csv_list=glob.glob('{}/nout*'.format(save_dir))
    num_csv_list = []
    if len(file_csv_list) != 0:
        for i in range(0, len(file_csv_list)):
            num_csv_list.append(int(file_csv_list[i].split('/')[-1].split('nout_')[1].split('.csv')[0]))
        max_num_csv_list.append(max(num_csv_list))
    else:
        max_num_csv_list.append(0)

    num_list = sorted(num_list)
    num_list_pop = [i for i in num_list if i <= min(max_num_csv_list)*every_num]
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
        if (number_out%every_num == 0):
            print('file_nout, csv_nout, numcsv = {}, {}, {}'.format(i, number_out, num_sta+i_data))
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
            df = df.assign(xcm = 1e+2*df['x_m'])
            data_nout_dict[i_data]=df
            i_data += 1

    return i_data, data_nout_dict, attr_nout_dict


if __name__ == "__main__":

    # print(sys.argv)
    
    ### setting params
    read_dir_name = './RESULTS_sph_flame_0.02mJ_hdf'

    every_num = 10
    every_num = 1

    read_path = '{}'.format(read_dir_name)
    file_list=glob.glob('{}/results*'.format(read_path))
    print(file_list)

    dir_name = read_dir_name.split('/')[-1].split('_hdf')[0]
    save_dir = dir_name+'_csv'
    print('save directory:', save_dir)

    num_list, filepath_dic = get_filenum()
    num_sta = int(num_list[0]/every_num)

    i_data = 0
    for i in num_list:
        if i == num_list[-1]:
            try:
                f = h5py.File(filepath_dic[i], 'r')
            except:
                continue
        i_data, data_nout_dict, attr_nout_dict = read_HDFfiles(i_data)

        for i_df, df in data_nout_dict.items():
            save_path = '{0}/t{1:1.6f}s_nout_{2:05}.csv'.format(save_dir, 
                        float(attr_nout_dict[i_df]['time from 0s']), num_sta+i_df)
            df.to_csv(save_path, index=False)


