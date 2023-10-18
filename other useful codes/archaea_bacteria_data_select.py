#!/usr/bin/env python
# coding: utf-8

import os
import re
from pyfaidx import Fasta
from tqdm import tqdm


# # 处理各古菌的fasta文件
# - 对前面挑选出来的古菌的fasta文件进行处理：
# - 挑选出每个fasta文件中每个氨基酸序列的头信息中‘description’部分为“hypothetical protein”的氨基酸序列；
# - 并将这些挑选出来的氨基酸单独存储成新的文件。


def parse_one_long_name(long_name:str):
    '''对1个氨基酸的长名称进行解析，这里主要判断长名称中，"description"部分是否为：hypothetical protein/hypothetical/possible predicted protein等'''
    
    marks = ['hypothetical protein','hypothetical','possible predicted protein']
    result = re.match('^.*description:(.*)$',long_name)  # 正则化匹配出长名称中的description内容，用于后续判断蛋白质类型
    try:
        # print(result.group(1).strip())  # 正则化匹配结果
        if result.group(1).strip() in marks:  # 判断描述信息中是否将该蛋白描述为“hypothetical protein”，是则将mark标记为True
            return 'Yes'
        else:
            return 'No'
    except:
        return 'Fail'


def pick_acids_accord_seq_names(proteins_obj, new_fasta_file, short_names_list, long_name=True, output_info=False):
    
    '''函数功能：从某个fasta文件中，根据氨基酸序列的“短名称”提取对应氨基酸序列，并存储成新的fasta文件；
        proteins_obj：Fasta(xxx.fasta)读取的fasta原文件对象；
        new_fasta_file：用于存储选中氨基酸的文件（存成fasta格式），命名时请以.fasta结尾；
        short_names_list：包含氨基酸序列名称的list, []；
        long_name：新的fasta文件中氨基酸信息是按长名称还是短名称存储。
    '''
    
    if not isinstance(short_names_list,list):
        raise ValueError("Warning:'short_names_list' should be a list type.")

    seq_fail = []  # 存放提取失败的氨基酸名称

    with open(new_fasta_file, 'w', encoding='utf-8') as f:

        for name_i in short_names_list:
            try:
                if long_name==True:
                    name = proteins_obj[name_i].long_name  # 从原始fasta文件中提取对应氨基酸短名称的长名称
                else:
                    name = proteins_obj[name_i].name  # 从原始fasta文件中提取对应氨基酸短名称的短名称（实际上就是name_i）
                # print(name)
                
                # 写入氨基酸名称
                f.write('>'+name.strip()+'\n')  

                # 写入氨基酸序列
                for seq_seg in proteins_obj[name_i]:  # 逐行提取某个蛋白质的氨基酸序列
                    f.write(str(seq_seg)+'\n')

            except:
                seq_fail.append(name_i)
    
    if output_info:
        print(f'New fiel:{new_fasta_file}; {len(short_names_list)-len(seq_fail)} proteins succeed, {len(seq_fail)} failed!')                  

    return seq_fail


# 主程序
if __name__=="__main__":
    # archaea_data_path = '/dload_bacteria_protein_data_cgx/archaea_data'  # 存储古菌数据文件的路径
    archaea_data_path = './archaea_data'  # 存储古菌数据文件的路径
    archaea_fasta_list = [file_i for file_i in os.listdir(archaea_data_path) if file_i[-3:]=='.fa']  # 该项目中fasta文件以.fa结尾（以防文件中还有其他格式文件）

    save_seleted_fasta_path = archaea_data_path+'/'+'selected_archaea_cgx'

    pick_seq_length = 400  # 选择的氨基酸长度

    try:
        os.makedirs(save_seleted_fasta_path)  # 创建路径，用于存储新建的fasta文件
    except:
        pass

    for archaea_fasta_i in tqdm(archaea_fasta_list):  # 不同fasta文件遍历
        archaea_fasta_path = archaea_data_path+'/'+archaea_fasta_i  # 包含路径信息的fasta文件名
        #print(archaea_fasta_path)
        proteins = Fasta(archaea_fasta_path)  # 读取fasta文件
        short_names = list(proteins.keys())  # 得到所有short names

        result = re.match('^(.*).fa',archaea_fasta_i)  # 截取原fasta文件的前面部分作为新文件名的一部分
        new_fasta_file = save_seleted_fasta_path+'/'+result.group(1)+'.cgxselect'+'.fasta'  # 加上“.cgxselect”作为新名称的标记

        pick_short_names = []  # 存储被选中氨基酸序列的short names
        for short_name_i in short_names:  # short name 循环
            long_name_i = proteins[short_name_i].long_name  # 根据短名称，获取长名称
            mark = parse_one_long_name(long_name_i)  # 对长名称进行判断，选中的返回True
            seq_len = proteins[short_name_i].unpadded_len  # 得到当前氨基酸序列的长度
            
            if (mark=='Yes') and (seq_len<=pick_seq_length):
                pick_short_names.append(short_name_i)  # 记录选中的氨基酸的短名称
            elif mark=='Fail':  # 出现异常
                print(f"Warning:Protein {short_name_i} in {archaea_fasta_i} get no description.")
            else:
                pass

        # 对当前fasta文件，根据选中氨基酸序列的短名称pick_short_names，从原序列中选择对应氨基酸序列并存储在新文件new_fasta_file中
        fail_amino = pick_acids_accord_seq_names(proteins, new_fasta_file, pick_short_names)

        # 删除因使用Fasta()方法而自动生成的.fai文件，对项目没啥用。
        os.remove(proteins.faidx.filename+'.fai')