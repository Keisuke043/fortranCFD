
import os
import sys
import glob
import subprocess
import numpy as np
import pandas as pd
import h5py


def detect_peaks(df, detect_param, extract_params):

    peak_dict = {key:[] for key in extract_params}
    if (df.loc[0, detect_param] > df.loc[1, detect_param]) \
        and (xlim[0] <= df.loc[0, 'xcm'] <= xlim[1]):
        for i_param, param_name in enumerate(extract_params):
            peak_dict[param_name].append(df.loc[0, param_name])
    for i in range(len(df.loc[:, detect_param])-2):
        if (df.loc[i, detect_param] < df.loc[i+1, detect_param] > df.loc[i+2, detect_param]) \
            and (xlim[0] < df.loc[i+1, 'xcm'] < xlim[1]):
            for i_param, param_name in enumerate(extract_params):
                peak_dict[param_name].append(df.loc[i+1, param_name])
    peak_list = []
    for param_name in extract_params:
        peak_list.append(peak_dict[param_name])

    return peak_list


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
            data_nout_dict[i_num_post]=df
            i_num_post += 1

    return i_num_post, data_nout_dict, attr_nout_dict


def write_peaks_csv():

    save_params = ['num','t_s']+extract_params
    num_last = filenum_list[-1]
    for i_peak in range(0, num_peaks):
        read_path_csv = read_path_csv_dic[i_peak]
        if not os.path.isfile(read_path_csv) or len(sys.argv) == 2:
            break
        f = open(read_path_csv, mode='r+')
        l_strip = [s.strip() for s in f.readlines()]
        f.close()
        l_last = l_strip[-1].split(',')
        df_last = pd.DataFrame(l_last, index=save_params).T
        num_last_tmp = int(df_last['num'])
        time_last_tmp = float(df_last['t_s'])
        if num_last_tmp < num_last:
            num_last = num_last_tmp

    if os.path.isfile(read_path_csv) and len(sys.argv) == 1:

        del_filenum = filenum_list[1] - filenum_list[0]
        write_filenum_list = [num for num in filenum_list if every_num*num_last - del_filenum < num]
        write_num_list = list(map(lambda filenum:int(filenum/every_num), write_filenum_list))

        for i_peak in range(0, num_peaks):
            read_path_csv = read_path_csv_dic[i_peak]
            df = pd.read_csv(read_path_csv)
            df = df[df['num']<write_num_list[0]]
            df.to_csv(read_path_csv, index=False)

    f_csv = []
    for i_peak in range(0, num_peaks):
        read_path_csv = read_path_csv_dic[i_peak]
        if os.path.isfile(read_path_csv):
            f = open(read_path_csv, mode='a')
        if not os.path.isfile(read_path_csv) or len(sys.argv) == 2:
            f = open(read_path_csv, mode='w')
            write_filenum_list = filenum_list
            write_num_list = list(map(lambda filenum:int(filenum/every_num), write_filenum_list))
            f.write('num,t_s')
            for param_name in extract_params:
                f.write(',{0}_{1}'.format(param_name, i_peak+1))
            f.write('\n')
        f_csv.append(f)

    num_sta = write_num_list[0]

    print(write_filenum_list)

    i_data = 0
    for i_numfile in write_filenum_list:
        if i_numfile == write_filenum_list[-1]:
            try:
                f = h5py.File(filepath_dic[i], 'r')
            except:
                continue
        i_data, data_nout_dict, attr_nout_dict = read_HDFfiles(i_numfile, i_data, num_sta)
        print('Writing in csv file {}'.format(i_data))
        for i_nout, df in data_nout_dict.items():
            list_peaks = detect_peaks(df, detect_param, extract_params)
            df_peaks = pd.DataFrame(list_peaks, index=extract_params).T
            df_peaks_sort = df_peaks.sort_values(detect_param, \
                            ascending=False).reset_index(drop=True)
            index_peaks = list(df_peaks_sort.index)
            for i_peak in range(0, num_peaks):
                time = float(attr_nout_dict[i_nout]['time from 0s'])
                f_csv[i_peak].write('{0},{1}'.format(num_sta+i_nout, time))
                if i_peak in index_peaks:
                    for variable_peak in df_peaks_sort.loc[i_peak]:
                        f_csv[i_peak].write(',{}'.format(variable_peak))
                    f_csv[i_peak].write('\n')
                else:
                    f_csv[i_peak].write('\n')
                f_csv[i_peak].flush()

        f_hdf.close()
    for i_peak in range(0, num_peaks):
        f_csv[i_peak].close()


if __name__ == "__main__":

    read_path_name = './RESULTS_sph_flame_*mJ_hdf'
    read_path_list = glob.glob(read_path_name)
    print(read_path_list)


    xlim = [0, 4.0] # cm
    detect_param = 'HRR'
    num_peaks = 3
    every_num = 10
    every_num = 1

    extract_params = 'all'
    extract_params = ['x_m', 'HRR']

    for read_path in read_path_list:

        file_list = glob.glob('{}/results_*.h5'.format(read_path))
        filenum_list, filepath_dic = get_filenum(file_list)
        f_hdf = h5py.File(filepath_dic[filenum_list[0]], 'r')

        print(read_path)
        print(filenum_list)

        save_dir = './RESULTS_sph_flame_peaks_ev{}_csv'.format(every_num)
        if extract_params == 'all':
            save_dir = save_dir + '_all'
            extract_params = list(f_hdf[list(f_hdf.keys())[0]].keys())

        ener_dir = read_path.split('_hdf')[0]
        save_dir = save_dir + '/{}'.format(ener_dir)

        read_path_csv_dic = {}
        for i_peak in range(0, num_peaks):
            read_path_csv = '{0}/peaks_{1}_{2}.csv'.format(save_dir, detect_param.strip(), i_peak+1)
            read_path_csv_dic[i_peak] = read_path_csv

        subprocess.call(['mkdir','-p', save_dir])
        write_peaks_csv()

        read_path = '{0}/peaks_{1}_*.csv'.format(save_dir, detect_param.strip())

        df_all = pd.DataFrame()
        for i_peak in range(0, num_peaks):
            read_path_csv = read_path_csv_dic[i_peak]
            df = pd.read_csv(read_path_csv)
            df_all = pd.concat([df_all, df], axis=1)
            df_all = df_all.loc[:,~df_all.columns.duplicated()]

        save_path = '{0}/peaks_{1}.csv'.format(save_dir, detect_param.strip())
        df_all.to_csv(save_path, index=False)


