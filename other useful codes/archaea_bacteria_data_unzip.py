#!/usr/bin/env python
# coding: utf-8

import os
import gzip
from tqdm import tqdm

# ## 解压古菌的压缩数据文件（.gz）
# - 解压出来是fasta格式文件

# 解压xx.gz压缩文件的函数
def un_gz(file_name):
     """ungz zip file"""
     f_name = file_name.replace(".gz", "")
     #获取文件的名称，去掉
     g_file = gzip.GzipFile(file_name)
     #创建gzip对象
     open(f_name, "wb+").write(g_file.read())
     #gzip对象用read()打开后，写入open()建立的文件里
     g_file.close()
     #关闭gzip对象


if __name__ == "__main__":
    archaea_data_path = '/dload_bacteria_protein_data_cgx/archaea_data'  # 存储古菌数据文件的路径  
    archaea_name_list = [file_i for file_i in os.listdir(archaea_data_path) if file_i[-3:]=='.gz']  # 得到压缩文件名（.gz结尾）
    for archaea_name_i in tqdm(archaea_name_list):
        archaea_name_path = archaea_data_path+'/'+archaea_name_i  # 包含路径信息的数据文件名
        try:
            un_gz(archaea_name_path)  # 解压xx.gz文件
        except:
            print(f'解压失败：{archaea_name_path}')