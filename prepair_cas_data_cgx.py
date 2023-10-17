import os
from rdkit import Chem
import time
import numpy as np
import random
import csv
from tqdm import tqdm
from matplotlib import pyplot as plt

if __name__ == '__main__':
    start = time.process_time()

###################################################################################################################################################
    most_caslen_2_smiles = 2000  # ★★★★★★★★★(经测试，在本机中，cas序列大于600时，用rdkit转化为smiles将导致内核奔溃，在服务器上不会)

    # 数据文件路径
    path = os.getcwd() + '/Cas_data'
    print(path)

    ###########################################################
    """
    将Cas蛋白质的氨基酸序列数据转化成smiles数据¶¶
    所有Cas蛋白的氨基酸序列数据以 *.fasta格式存储，一个文件对应1个Cas蛋白
    程序将每个Cas氨基酸序列转换成了对应的smiles格式数据
    """
    cas_len_tongji = np.zeros([10])  # 用于统计cas蛋白的序列长度范围
    len_one_bar = 200  # 每200个序列长度作为一个分割点，200、400、600……
    cas_most_length = 0  # 存放最长cas的序列长度

    num_cas2smi = 0  # 成功被转化成smiles序列的cas蛋白数

    smiles_len_tongji = []  # 存放每一个smiles序列的长度

    for cas_file in tqdm(os.listdir(path)):
        #     print(cas_file+'系列', flush=True)

        cas_file_path = path + '/' + cas_file  # 当前cas系列文件路径
        cas_data_names = os.listdir(cas_file_path)  # 当前cas系列文件路径下所有文件名（包括cas数据文件和其他文件）
        cas_data_names = [i for i in cas_data_names if i[-6:] == '.fasta']  # 筛选出cas数据文件，以“.fasta”结尾的文件

        cas_smiles = cas_file_path + '/' + cas_file + '_' + 'smiles' + '.csv'  # 在某一系列cas路径下建立一个csv文件，用于存储每个cas对应的的smiles序列
        with open(cas_smiles, "w", newline="") as datacsv:  # 新建csv文件
            # dialect为打开csv文件的方式，默认是excel，delimiter="\t"参数指写入的时候的分隔符
            csvwriter = csv.writer(datacsv, dialect=("excel"))
            # #写入第一行信息，即列名称。
            # csvwriter.writerow(["smiles","if_cas"])

            # 开始对某个cas系列下的每一个cas蛋白对应的fasta文件转化为对应的smiles序列，并写进同一个csv文件
            for data_i in range(len(cas_data_names)):
                cas_data_path = cas_file_path + '/' + cas_data_names[data_i]  # 带绝对路径的具体cas数据名
                # print(cas_data_path)
                with open(cas_data_path, 'r') as f:
                    cas_seq = ''
                    for line in f.readlines()[1:]:  # fasta数据的第1行为解释信息
                        # print(line.strip(), flush=True)
                        cas_seq += line.strip()  # 去除换行符并拼接
                    # print(cas_seq)
                    cas_len_tongji[int(np.ceil(len(cas_seq) / len_one_bar))] += 1  # cas蛋白序列长度位于哪个区段，对应区段加1（200作为1个间隔）

                    # 找出最长cas序列长度
                    if len(cas_seq) > cas_most_length:
                        cas_most_length = len(cas_seq)

                    if 0 < len(cas_seq) < most_caslen_2_smiles:  # ★★如果序列大于600将出现内核奔溃★★★★★★★★★★★★★★
                        mol = Chem.MolFromFASTA(cas_seq)  # 将序列转化成smiles格式数据
                        if mol is not None:
                            # print(f'转换开始！{cas_data_path}，序列长度：{len(cas_seq)}', flush=True)
                            cas_smiles = Chem.MolToSmiles(mol)
                            num_cas2smi += 1  # 记录转换成功的
                            smiles_len_tongji.append(len(cas_smiles))  # 记录当前smiles序列长度
                            csvwriter.writerow([cas_smiles, '1'])  # 将smiles格式数据保存到前面创建的csv文件中,并同时指定标签1

    print(f'总共{int(cas_len_tongji.sum())}个cas蛋白，转换成功{num_cas2smi}个！')
    print(f'cas最长长度：{cas_most_length}')
    print(f'smiles最长长度：{max(smiles_len_tongji)}')

    ###绘制统计图
    x_ticks = np.zeros([len(cas_len_tongji)])
    for i in range(len(cas_len_tongji)):
        x_ticks[i] = (i + 1) * len_one_bar
        print(f'cas length {x_ticks[i]}: {int(cas_len_tongji[i])}个cas')

    # 绘图
    fig, ax = plt.subplots(1, 1)
    plt.plot(x_ticks, cas_len_tongji)
    ax.set_xlabel('cas length', fontdict=dict(fontsize=16, color='r'))
    ax.set_ylabel('number of cas', fontdict=dict(fontsize=16, color='r'))
    plt.tick_params(axis='x', colors='red')
    plt.tick_params(axis='y', colors='red')
    plt.show()

    plt.hist(smiles_len_tongji, bins=20)
    plt.xlabel('smiles length', fontdict=dict(fontsize=16, color='r'))
    plt.ylabel('number of smiles', fontdict=dict(fontsize=16, color='r'))
    plt.tick_params(axis='x', colors='red')
    plt.tick_params(axis='y', colors='red')

    plt.show()

# ###################################################################################################################################################
#     """正负Cas蛋白数据合并形成训练集¶"""
#
#     # 读取非Cas蛋白smiles数据
#     nocas_smiles_path = 'Cas_data_nega/total_nocas_smiles_.csv'  # 非cas蛋白smiles存储路径
#     with open(nocas_smiles_path, encoding='utf-8') as nocasf:
#         nocas_smiles_data = nocasf.readlines()
#     count = len(nocas_smiles_data)  # 计算总共有多少个no cas蛋白
#
#     # 每一个类型的cas蛋白路径循环（Cas1、Cas2、……）
#     for cas_file in tqdm(os.listdir(path)):
#         # print(cas_file+'系列', flush=True)
#
#         cas_smiles_path = path + '/' + cas_file + '/' + cas_file + '_smiles.csv'  # cas数据文件路径
#         # print(cas_smiles_path)
#
#         # 待生成新文件（路径+名称），用于存储cas和nocas混合组成的smiles格式文件，作为有正负样本的样本集
#         cas_nocas_smiles_path = path + '/' + cas_file + '/' + cas_file + '_with_nocas_smiles.csv'
#         # print(cas_nocas_smiles_path)
#         with open(cas_nocas_smiles_path, "w", newline="") as datacsv:  # 生成对应的csv文件
#             # dialect为打开csv文件的方式，默认是excel，delimiter="\t"参数指写入的时候的分隔符
#             csvwriter = csv.writer(datacsv, dialect=("excel"))
#             # 写入第一行信息，即列名称。
#             csvwriter.writerow(["smiles", "if_cas"])
#
#             # 打开已经存在的cas蛋白的smiles文件（.csv）
#             with open(cas_smiles_path, encoding='utf-8') as f:
#                 for smiles_i in f.readlines():
#                     smiles_i = smiles_i.strip()
#                     # print(smiles_i[:-2])  # 抽出smiles格式数据
#                     # print(smiles_i[-1])  # 抽出1位标志0
#                     csvwriter.writerow([smiles_i[:-2], '1'])  # 将cas的smiles数据及标签“1”写入新的文件
#
#                     # 在nocas的smiles文件中随机提取一个smiles文件作为负样本
#                     nocas_smiles_i = nocas_smiles_data[np.random.randint(1, count)][:-3]
#                     # print(nocas_smiles_i)
#                     csvwriter.writerow([nocas_smiles_i, '0'])  # 将非cas的smiles数据及标签“0”写入新的文件



###################################################################################################################################################

    end = time.process_time()
    runTime = end - start
    print("运行时间：", runTime)