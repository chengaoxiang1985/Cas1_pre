## 从蛋白质下载链接上下载蛋白质数据（在前述保存信息的技术上独立运行）
import os
import requests
import re
import json
import pandas as pd
import multiprocessing
import time

# 读取前面保存的存储的所有细菌的蛋白质相关信息文件（csv文件）
file_name = 'df_protein_info.csv'
df = pd.read_csv(file_name, encoding='utf-8')

# 先获取A-Z和other路径中已经下载好的数据，避免重复下载浪费时间，存储在字典data_file_dict中
data_path_l1 = './download_bacteria_protein_data'  # 存放A-Z和other目录的目录
data_path_l2 = os.listdir(data_path_l1)  # 获取A-Z和other目录
data_names = [os.listdir(data_path_l1+"/"+data_file_name) for data_file_name in data_path_l2]  # 获取A-Z和other下已经存在的实际数据文件名称
data_file_dict = dict(zip(data_path_l2, data_names))
# print(data_file_dict.keys())
# print(data_file_dict)

def Download_one_data(bacteria_name, data_file_name, data_url, save_path):
    '''下载一个压缩蛋白数据
        bacteria_name：细菌名称
        data_file_name：delattr白压缩数据名称
        data_url：蛋白压缩数据url
        save_path：数据保存路径
    '''
    
    fail_bacteria = []
    try:
        rest = requests.get(data_url)
        with open(save_path+'/'+data_file_name,'wb') as file:
            file.write(rest.content)
    except:
        fail_bacteria.append(bacteria_name)
        print
    
    return fail_bacteria

# 用于多线程下载 
def Download_multiprocess(i):

    first_cal = df['细菌名称'][i][0]  # 获取细菌名称的首字符（可能是大写字母、小写字母、其他符号）
    
    if not first_cal.isalpha():  # 如果首字符不是字母，则数据存储在other文件中
        data_path_l2_name = 'other' 
    elif first_cal in 'abcdefghijklmnopqrstuvwxyz':  # 部分细菌名称的首字符是小写字母，都将其存储在大写字母的路径中（A-Z）
        first_cal = first_cal.upper()  # 小写字母换成大写，便于在对应大写字母文件中查找
        data_path_l2_name = first_cal
    else:  # 其他情况都是大写字母开头
        data_path_l2_name = first_cal
    
    
    # 得到对应首字母的文件夹
    if first_cal in list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'): 
        save_path = './download_bacteria_protein_data'+'/'+first_cal
    else:
        save_path = './download_bacteria_protein_data'+'/'+'other'
    # print(save_path)

    if df['蛋白质压缩数据名称'][i] not in data_file_dict[data_path_l2_name]:  # 判断当前待下载的蛋白质数据是否已经存在了，如果存在则无需下载

        fail_bacteria = Download_one_data(df['细菌名称'][i], df['蛋白质压缩数据名称'][i], df['蛋白质下载页'][i], save_path)

        if fail_bacteria:
            fail_bacteria_names.appenda(fail_bacteria)

        print(f'共{df.shape[0]}个，已下载：{i+1}个！，存储路径：{save_path}')
    else:
        current_data_name = df['蛋白质压缩数据名称'][i]
        print(f'已存在:{current_data_name}')


if __name__=='__main__':
    
    print('Program starting......')
    print(time.strftime( '%Y-%m-%d %H:%M:%S' ,time.localtime(time.time( ))))
    start_time = time.time( )
    
    fail_bacteria_names = []  # 存放没有被成功下载的数据名称
    
    # # choice 1 挨个下载，慢###################################################################################################
    # #  挨个下载
    # for i in range(df.shape[0]):  # 所有细菌循环
    #     first_cal = df['细菌名称'][i][0]  # 获取细菌名称的首字母
    #     # 得到对应首字母的文件夹
    #     if first_cal in list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'): 
    #         save_path = './download_bacteria_protein_data'+'/'+first_cal
    #     else:
    #         save_path = './download_bacteria_protein_data'+'/'+'other'
    #     # print(save_path)

    #     fail_bacteria = Download_one_data(df['细菌名称'][i], df['蛋白质压缩数据名称'][i], df['蛋白质下载页'][i], save_path)

    #     if fail_bacteria:
    #         fail_bacteria_names.append(fail_bacteria)

    #     print(f'共{df.shape[0]}个，已下载：{i+1}个！，存储路径：{save_path}')


    # choice 2 多进程下载，快###################################################################################################
    pool = multiprocessing.Pool(8)  # 开启进程数
    pool.map(Download_multiprocess, [i for i in range(df.shape[0])])
    pool.close()  # 关闭进程池，不再接受新的进程
    pool.join()  # 主进程阻塞等待子进程的退出


    end_time = time.time( )
    print(f'Total running time : {end_time - start_time} s')