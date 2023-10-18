# ## 说明：从原始fasta文件中，根据数量选取氨基酸序列，并存储成新的fasta文件。（包括随机选取、等间隔选取）

from typing import List
from pyfaidx import Fasta
import random
import numpy as np
import os
from matplotlib import pyplot as plt


def pick_acids_accord_seq_number(raw_fasta_file, new_fasta_file, acids_select_num, choice='random', random_seed=2, long_name=True):
    '''Get acids seqs from fasta file according position index and saved as a new fasta file.
        raw_fasta_file: raw fasta file, where to select acids seqs
        new_fasta_file: new fasta file to save the selected acids seqs
        acids_select_num: int, how many acids to be selected
        choice: ='random' random select (no duplicate ); ='equal' equally spaced selection
        random_seed: set seed for random.sample()
    '''

    if choice not in ['random','equal']:
        raise NameError('Parameter "choice" should be "random" or "equal".')

    proteins_obj = Fasta(raw_fasta_file)  # 读取fasta文件
    os.remove(proteins_obj.filename+'.fai')  # 删除因使用Fasta()方法而自动生成的.fai文件，对项目没啥用。

    acids_names = list(proteins_obj.keys())
    acids_num = len(acids_names)

    if acids_select_num>=acids_num:  # if acids_select_num is larger than the total acids number of the raw fasta file, then no need to select.  
        print(f'Warning: acids_select_num={acids_select_num} is larger than the total acids number={acids_num} in the raw fasta file!')
        print('So there is no need to select!')
    else:
        if choice=='random':
            random.seed(random_seed)
            # randomly and unrepeatably select 'acids_select_num' numbers within [0,acids_num), and used as acids position index to be selected
            position_index_list = random.sample(range(acids_num),acids_select_num)  
        elif choice=='equal':
            inter = int(np.floor(acids_num/acids_select_num))  # Position interval
            position_index_list = [i*inter for i in range(acids_select_num)]

        position_index_list = list(set(position_index_list))  # find the duplicated position index and remove

        failed_index = []

        with open(new_fasta_file, 'w', encoding='utf-8') as f:

            for p_i in position_index_list: # loop of position index 
                try:
                    acid_name_i = acids_names[p_i]  # get the short name of the acid at p_i

                    if long_name is True:
                        name = proteins_obj[acid_name_i].long_name  # 从原始fasta文件中提取对应氨基酸短名称的长名称
                    else:
                        name = proteins_obj[acid_name_i].name

                    seq = proteins_obj[acid_name_i][:].seq  # get the acid sequence 
                    f.write('>'+name.strip()+'\n')  # 写入氨基酸名称
                    one_seq_frgs = segment_one_seq(seq)  # 将当前氨基酸分割成等长片段
                    for frg in one_seq_frgs:
                        f.write(frg+'\n')  # 写入氨基酸

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


    new_fasta_file = 'select_accord_seq_number.FASTA'

    failed_index = pick_acids_accord_seq_number(raw_fasta_file, new_fasta_file, acids_select_num=4, choice='random', random_seed=2, long_name=True)  # 核心计算，从原fasta文件中提取氨基酸序列
    print(f"The following position index are failed to be selected:{failed_index}")

    new_proteins = Fasta(new_fasta_file)  # 读取新生成的fasta文件
    os.remove(new_proteins.filename+'.fai')  # 删除因使用Fasta()方法而自动生成的.fai文件，对项目没啥用。

    leng_list = []
    for name_i in new_proteins.keys():
        leng_list.append(len(new_proteins[name_i][:].seq))
    show_len_distribution_by_histgram(leng_list, new_fasta_file+' sequence length distribution')