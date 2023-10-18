# # 说明：该函数从1个原始的fasta文件中，根据氨基酸的长度挑选出满足长度要求的氨基酸，并存储成新的fasta格式文件。

from pyfaidx import Fasta
from matplotlib import pyplot as plt
import os


# # 定义函数
def pick_acids_accord_seq_lenth(raw_fasta_file, new_fasta_file, min_seq_len, max_seq_len, long_name=True):
    '''
        raw_fasta_file: raw fasta file, where to select acids seqs
        new_fasta_file：用于存储选中氨基酸的文件（存成fasta格式），命名时应以.fasta结尾
        min_seq_len：要求的最小氨基酸长度
        max_seq_len：要求的最大氨基酸长度
        long_name: 是否存储长名字
    '''
    
    proteins_obj = Fasta(raw_fasta_file)
    os.remove(proteins_obj.filename+'.fai')  # 删除对应fai文件。
    
    seq_succeed = 0
    with open(new_fasta_file, 'w', encoding='utf-8') as f:
        
        for name_i in proteins_obj.keys():
            try:
                if long_name==True:
                    name = proteins_obj[name_i].long_name  # 从原始fasta文件中提取对应氨基酸短名称的长名称
                else:
                    name = proteins_obj[name_i].name

                if min_seq_len<=proteins_obj[name_i].unpadded_len<=max_seq_len:  # 判断当前氨基酸长度是否满足要求
                    seq_succeed += 1  # 统计被选中的氨基酸个数
                    f.write('>'+name.strip()+'\n')  # 写入氨基酸名称
                    
                    for seq_seg in proteins_obj[name_i]:  # 逐行提取某个蛋白质的氨基酸序列
                        f.write(str(seq_seg)+'\n')  # 写入氨基酸
                        
            except:
                print(f'{name_i} Failed......')

    print(f'New file has been created! {seq_succeed} are selected!')

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


# # 使用

if __name__ == '__main__':
    raw_fasta_file = 'test.fasta'

    proteins = Fasta(raw_fasta_file)  # 读取fasta文件
    os.remove(proteins.filename+'.fai')  # 删除对应fai文件。

    leng_list = []
    for name_i in proteins.keys():
        leng_list.append(len(proteins[name_i][:].seq))
    show_len_distribution_by_histgram(leng_list, 'All sequence length distribution')


    # 根据氨基酸长度挑选，并存成新文件
    new_fasta_file = 'select_accord_seq_length.FASTA'  # 新fasta文件名称

    pick_acids_accord_seq_lenth(raw_fasta_file, new_fasta_file, min_seq_len=200, max_seq_len=400, long_name=True)  # 根据序列选择函数

    new_proteins = Fasta(new_fasta_file)  # 读取新生成的fasta文件
    os.remove(new_proteins.filename+'.fai')  # 删除对应fai文件。

    leng_list = []
    for name_i in new_proteins.keys():
        leng_list.append(len(new_proteins[name_i][:].seq))
    show_len_distribution_by_histgram(leng_list, new_fasta_file+' sequence length distribution')