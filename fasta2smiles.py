#!/usr/bin/env python
# coding: utf-8

# # 将fasta文件转化成smiles格式文件

from rdkit import Chem
import csv
from tqdm import tqdm
from pyfaidx import Fasta
import os


# # 定义函数
def Fasta2Smiles(fasta_name, smiles_name, label, seq_len_most, smiles_len_most=10000000):
    '''amino acid sequence in fasta to smiles
       fasta_name: fasta file name
       smiles_name: smiles file nat
       label: point if Cas seq or not， 1 for Cas, 0 for non Cas
       seq_len_most: sequence with length more than this number will not be transform to smiles
       smiles_len_most: smiles with length less than this number will be stored into new file
    '''
    
    total_acid = 0  # 总共过少氨基酸序列
    num_seq2smi = 0  # 记录转换成功的氨基酸
    smiles_len = []  # 记录每个smile
    acids_names_select = []  # 存储被选中氨基酸的名称
    
    with open(smiles_name, "w", newline="") as datacsv:  # 新建csv文件
        csvwriter = csv.writer(datacsv, dialect=("excel"))
        #写入第一行信息，即列名称。
        csvwriter.writerow(["smiles","if_cas1"])

        proteins = Fasta(fasta_name)  # 打开原始蛋白文件
        os.remove(proteins.filename+'.fai')  # 删除对应fai文件。
        names = list(proteins.keys())  # 获取所有蛋白名称
        for name_i in tqdm(names):  # Fasta自动读取fasta文件并处理
            total_acid += 1
            seq = proteins[name_i][:].seq
            if 0 < len(seq) <= seq_len_most: 
                mol = Chem.MolFromFASTA(seq)  # 将序列转化成smiles格式数据
                if mol is not None:
                    seq_smiles = Chem.MolToSmiles(mol)
                    if len(seq_smiles)<=smiles_len_most:
                        num_seq2smi += 1  # 记录转换成功的
                        smiles_len.append(len(seq_smiles))  # 记录当前smiles序列长度
                        csvwriter.writerow([seq_smiles, label])  # 将smiles格式数据保存到前面创建的csv文件中,并同时指定标签1
                        acids_names_select.append(name_i)  # 将选中的氨基酸名称记录下来
    print('Smiles file is done!')
    
    return total_acid, num_seq2smi, smiles_len, acids_names_select


# # 使用
if __name__== "__main__":

    # 综合参数：
    seq_len_most=30000  # ★★在PC上运行时，如果氨基酸序列大于600将出现内核奔溃（服务器上没问题）★，在服务器上OK
    smiles_len_most=80000  # 选择smiles格式长度小于8000的氨基酸

    # # 转换数据文件1###############################################################################################################
    # fasta_name1 = 'uniref50_len400to1300.fasta'
    # smiles_name1 = 'uniref50_len400to1300_smiles.csv'
    # Cas1_total_acid, Cas1_num_seq2smi, Cas1_smiles_len, acids_names_select =  Fasta2Smiles(fasta_name1, 
    #                                                                                         smiles_name1, 
    #                                                                                         label=0,
    #                                                                                         seq_len_most=seq_len_most,
    #                                                                                         smiles_len_most=smiles_len_most)

    # print(f'共{Cas1_total_acid}个Cas1，成功转换{Cas1_num_seq2smi}个！')
    # print(f'最长smiles：{max(Cas1_smiles_len)}')


    # # 转换数据文件2###############################################################################################################
    # fasta_name1 = 'Cas1_len400to1300.FASTA'
    # smiles_name1 = 'Cas1_len400to1300_smiles.csv'
    # Cas1_total_acid, Cas1_num_seq2smi, Cas1_smiles_len, acids_names_select =  Fasta2Smiles(fasta_name1, 
    #                                                                                         smiles_name1, 
    #                                                                                         label=1,
    #                                                                                         seq_len_most=seq_len_most,
    #                                                                                         smiles_len_most=smiles_len_most)
    # print(f'共{Cas1_total_acid}个Cas1，成功转换{Cas1_num_seq2smi}个！')
    # print(f'最长smiles：{max(Cas1_smiles_len)}')


    # 转换数据文件3###############################################################################################################
    fasta_name1 = 'Cas1_combine.FASTA'
    smiles_name1 = 'Cas1_combine_smiles.csv'
    Cas1_total_acid, Cas1_num_seq2smi, Cas1_smiles_len, acids_names_select =  Fasta2Smiles(fasta_name1, 
                                                                                            smiles_name1, 
                                                                                            label=1,
                                                                                            seq_len_most=seq_len_most,
                                                                                            smiles_len_most=smiles_len_most)
    print(f'共{Cas1_total_acid}个Cas1，成功转换{Cas1_num_seq2smi}个！')
    print(f'最长smiles：{max(Cas1_smiles_len)}')