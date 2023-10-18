# ## 说明：从原始fasta文件中，根据序号（用户指定）切取对应位置的氨基酸序列，并存储成新的fasta文件。

from typing import List
from pyfaidx import Fasta
import os
from matplotlib import pyplot as plt


# ## 函数
def get_acids_accord_index(raw_fasta_file, new_fasta_file, position_index_list, long_name=True):
    '''Get acids seqs from fasta file according position index and saved as a new fasta file.
        raw_fasta_file: raw fasta file, where to select acids seqs
        new_fasta_file: new fasta file to save the selected acids seqs
        position_index_list: an int list that gives the position index of the target acids 
    '''
    
    proteins = Fasta(raw_fasta_file)  # 读取fasta文件
    os.remove(proteins.filename+'.fai')  # 删除因使用Fasta()方法而自动生成的.fai文件，对项目没啥用。
    
    acids_names = list(proteins.keys())
    acids_num = len(acids_names)

    position_index_list = list(set(position_index_list))  # find the duplicated position index and remove

    failed_index =[]

    with open(new_fasta_file, 'w', encoding='utf-8') as f:

        for p_i in position_index_list: # loop of position index 
            try:
                acid_name_i = acids_names[p_i]  # get the short name of the acid at p_i

                if long_name==True:
                    name = proteins[acid_name_i].long_name  # 从原始fasta文件中提取对应氨基酸短名称的长名称
                else:
                    name = proteins[acid_name_i].name

                # 写入氨基酸名称
                f.write('>'+name.strip()+'\n')  # 写入氨基酸名称
                # 写入氨基酸序列
                for seq_seg in proteins[acid_name_i]:  # 逐行提取某个蛋白质的氨基酸序列
                    f.write(str(seq_seg)+'\n')

            except:
                failed_index.append(p_i)

    print('New fasta file has been created!')
    if failed_index:
        print(f'{len(failed_index)} faild!')

    return failed_index


def show_len_distribution_by_histgram(list_num, fig_title):
    '''show squence length distribution by histgram
       list_num: a list of numbers or a one-dimensional numpy array.
    '''

    plt.figure(figsize=(30, 8))
    # plt.hist(list_num,range=(0,max(list_num)),bins=200,edgecolor='black',log=True)  # 带log方式显示
    plt.hist(list_num,range=(0, max(list_num)),bins=200,edgecolor='black')
    plt.tick_params(axis='x' ,colors= 'red' ,labelsize=20)
    plt.tick_params(axis='y' ,colors= 'red' ,labelsize=20)
    plt.xlabel('sequence length', fontdict=dict(fontsize=20, color='r'))
    plt.ylabel('number of sequence', fontdict=dict(fontsize=20, color='r'))
    plt.title(fig_title, fontdict=dict(fontsize=20, color='r'))
    plt.show()


# ## 使用（测试）
if __name__ == '__main__':
    raw_fasta_file = 'test.fasta'  # 测试数据

    proteins = Fasta(raw_fasta_file)  # 读取fasta文件
    os.remove(proteins.filename+'.fai')  # 删除因使用Fasta()方法而自动生成的.fai文件，对项目没啥用。

    leng_list = []
    for name_i in proteins.keys():
        leng_list.append(len(proteins[name_i][:].seq))
    show_len_distribution_by_histgram(leng_list, 'All sequence length distribution')

    new_fasta_file = 'select_accord_seq_index.FASTA'

    position_index_list = [2,1,2,3,4,'a']  # 实际中，位置序列需要自己生成★★★★

    failed_index = get_acids_accord_index(raw_fasta_file, new_fasta_file, position_index_list, long_name=True)  # 核心计算，从原fasta文件中提取氨基酸序列
    print(failed_index)


    new_proteins = Fasta(new_fasta_file)  # 读取新生成的fasta文件
    os.remove(new_proteins.filename+'.fai')  # 删除因使用Fasta()方法而自动生成的.fai文件，对项目没啥用。

    leng_list = []
    for name_i in new_proteins.keys():
        leng_list.append(len(new_proteins[name_i][:].seq))
    show_len_distribution_by_histgram(leng_list, new_fasta_file+' sequence length distribution')