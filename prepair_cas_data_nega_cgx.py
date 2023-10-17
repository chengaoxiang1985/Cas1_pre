import os
from rdkit import Chem
import time
import numpy as np
import csv
from tqdm import tqdm

"""
将非Cas蛋白质的氨基酸序列数据转化成smiles数据
# 原始数据由陈觅拷贝而来，氨基酸序列数据以*.faa格式存储
# 程序将每个氨基酸序列转换成了对应的smiles格式数据
# 注意：在pc上当氨基酸长度大于600后因内存问题而奔溃，在服务器上不会
"""

# 数据文件路径
path = os.getcwd()+'/Cas_data_nega'
print(path)

def data_file_to_list(data_file_name):
    """将一个包含多个氨基酸序列的文件进行处理
    data_file_name：数据文文件名
    return results_list：包含所有氨基酸序列的list（一个序列对应一个蛋白质）
    """

    results_list = []  # 存储所有蛋白的完整氨基酸序列
    with open(data_file_name, encoding='utf-8') as f:

        protein_seq = ''  #存储单个氨基酸序列
        for line in f.readlines():
            if line[0] == '>':  # 判断当前行是否为氨基酸序列的说明信息行
                if protein_seq != '':
                    results_list.append(protein_seq)  # 将单个蛋白的完整氨基酸序列存入results_list
                protein_seq = ''  # 当前氨基酸序列清零，用于存储下一个氨基酸序列
                continue  #跳出当前行
            else:  # 如果当前行是氨基酸序列行
                line = line.strip() # 去除当前行的换行符
                line = line[:-1] if line[-1] == "*" else line #去除序列尾部*

                protein_seq += line  # 同一个蛋白质的氨基酸序列进行拼接
    return results_list


if __name__ == '__main__':
    start = time.process_time()

    ######################################################################################################################################################
    # # 将非Cas蛋白质的氨基酸序列数据转化成smiles数据
    # num_nocas2smi = 0  # 成功被转化成smiles序列的cas蛋白数
    #
    # most_nocaslen_2_smiles = 2000  # ★★★★★★★★★(经测试，在本机中，cas序列大于600时，用rdkit转化为smiles将导致内核奔溃)
    #
    # for nocas_data_name_i in tqdm(os.listdir(path)):
    #
    #     if nocas_data_name_i[-3:] == 'faa':  # 如果是数据文件，末尾字符通常是‘faa’
    #         print(f'当前数据名：{nocas_data_name_i}', flush=True)  # 单纯的数据文件名
    #         nocas_smiles_file_i = path + '/' + nocas_data_name_i + '_smiles' + '.csv'  # 转化成csv文件后的全路径
    #         print(f'待生成csv文件名：{nocas_smiles_file_i}', flush=True)
    #
    #         with open(nocas_smiles_file_i, "w", newline="") as datacsv:  # 生成对应的csv文件
    #             # dialect为打开csv文件的方式，默认是excel，delimiter="\t"参数指写入的时候的分隔符
    #             csvwriter = csv.writer(datacsv, dialect=("excel"))
    #
    #             current_seq_list = data_file_to_list(path + '/' + nocas_data_name_i)  # 打开当前data文件并返回由每个蛋白的氨基酸序列组成的list
    #             # print(current_seq_list)
    #             for nocas_seq in current_seq_list:
    #                 if len(nocas_seq) < most_nocaslen_2_smiles:  # ★★如果序列大于600将出现内核奔溃（服务器上不会）★★★★★★★★★★★★★★
    #                     mol = Chem.MolFromFASTA(nocas_seq)
    #                     if mol is not None:
    #                         # print(f'转换开始！{cas_data_path}，序列长度：{len(cas_seq)}', flush=True)
    #                         nocas_smiles = Chem.MolToSmiles(mol)
    #                         num_nocas2smi += 1  # 记录转换成功的
    #                         csvwriter.writerow([nocas_smiles, "0"])  # 将smiles格式数据保存到前面创建的csv文件中
    #                         # print(cas_smiles, flush=True)


######################################################################################################################################################
    """对分割的csv文件进行合并处理;
    1个csv文件存储了部分smiles文件;
    将所有分割的非Cas蛋白的smiles文件;
    """
    nocas_smiles = [file for file in os.listdir(path) if file[-11:] == '_smiles.csv']  # 找出当前文件夹下所有以“_smiles.csv”结尾的文件

    if nocas_smiles is not []:  # 如果非空（即文件存在）
        count = 0

        file_name = 'total_nocas_smiles_.csv'  # 待生成的综合csv文件名称
        with open(path + '/' + file_name, "w", newline="") as datacsv:  # 生成对应的csv文件
            # dialect为打开csv文件的方式，默认是excel，delimiter="\t"参数指写入的时候的分隔符
            csvwriter = csv.writer(datacsv, dialect=("excel"))

            for file_i in tqdm(nocas_smiles):  # 逐个打开csv文件
                data_path = path + "/" + file_i
                print(f'当前数据名：{data_path}', flush=True)  # 单纯的数据文件名

                with open(data_path, encoding='utf-8') as f:
                    for smiles_i in f.readlines():
                        count += 1
                        smiles_i = smiles_i.strip()
                        # print(smiles_i[:-2])  # 抽出smiles格式数据
                        # print(smiles_i[-1])  # 抽出1位标志0
                        csvwriter.writerow([smiles_i[:-2], smiles_i[-1]])

        print(f'共处理{count}条smiles数据！', flush=True)

    else:
        print('没找到以“_smiles.csv”结尾的文件！')
    ######################################################################################################################################################



    end = time.process_time()
    runTime = end - start
    print("运行时间：", runTime)