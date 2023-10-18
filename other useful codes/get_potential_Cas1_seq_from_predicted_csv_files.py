import os
import pandas as pd
import re
from pyfaidx import Fasta
from tqdm import tqdm

def Get_potential_Cas_from_prediction(predicted_csv, score_start:float, score_end:float):
    """从1个预测结果文件（csv文件）中挑选出预测结果分数介于score_start和score_end之间的蛋白信息。
       返回被选中的蛋白在fasta文件中的名称，以及对应的预测值。
    """
    if (score_start<0 or score_start>1) or (score_end<0 or score_end>1):
        raise ValueError('score_start and score_end shoul be between 0 and 1!')
    if score_end<score_start:
        raise ValueError('score_end should equal or larger than score_start!')

    df = pd.read_csv(predicted_csv)  # 读取已经存在的数据文件
    pd_select = df[(score_start<=df['if_cas']) & (df['if_cas']<=score_end)]  # 筛选数据

    return pd_select


if __name__=="__main__":

    print('Program starting......')

    # 主要参数设置
    score_start=0.9  # 选择预测分数大于多少的蛋白质作为潜在的cas1蛋白
    score_end=1.0  # 选择预测分数大于多少的蛋白质作为潜在的cas1蛋白

    potential_cas1_fasta_file = 'POTENTIAL_CAS1_score_'+str(score_start)+'_'+str(score_end)+'.FASTA'  # 最后生成的fasta文件名
    # print(potential_cas1_fasta_file)
    potential_cas1_fasta_file_path = os.path.join(os.getcwd(),potential_cas1_fasta_file)     # 最后生成的fasta文件的存放路径
    # print(potential_cas1_fasta_file_path)

    predict_result_path = './selected_archaea_cgx'  # 存放预测结果文件的路径
    # os.listdir(predict_result_path)
    
    predict_file_name = [i for i in os.listdir(predict_result_path) if i[-12:]=='.predict.csv']  # 提取预测结果文件（csv文件）名称list

    ##################################################################################################################################
    # 将被选中的潜在cas1所在的fasta文件的名称,对应潜在cas1名和预测值，从各个csv文件中挑选出，并存储在字典中，形如：
    # {fasta1：[[cas1,cas1 cas1…],[sco,sco,sco…]]，fasta2：[[cas1,cas1，cas…],[sco,sco,sco…]]，……}
    selected_Cas_fasta_dict = {}
    for csv_i in tqdm(predict_file_name):  # 多个预测的csv文件遍历
        predicted_csv_path = os.path.join(predict_result_path, csv_i)

        # 根据csv文件名，生成对应的fata文件名（最终选中的cas蛋白的原始信息要到fasta文件中去取）
        match_res = re.match('^(.*).predict.csv',csv_i)
        fasta_file_name = match_res.group(1)+'.fasta'

        # 从一个预测的csv文件中挑选可能得cas
        pd_select = Get_potential_Cas_from_prediction(predicted_csv_path, score_start=score_start, score_end=score_end)

        # 将fasta文件名，选中的cas和对应的score作为字典的一个元素插入
        if not pd_select.empty:  # 只处理非空结果（即找到的）
            selected_Cas_fasta_dict[fasta_file_name] = [list(pd_select['protein_name']),list(pd_select['if_cas'])]


    ##################################################################################################################################
    # 根据前面生成的字典，从原始fasta文件中奖对应潜在cas1的信息和序列抽出，并存储在一个新的fasta文件中。

    # potential_cas1_fasta_file_path = os.path.join(predict_result_path, potential_cas1_fasta_file)
    # print(potential_cas1_fasta_file_path)

    num_of_successes = 0  # 记录成功了多少个氨基酸序列

    f = open(potential_cas1_fasta_file_path, 'w', encoding='utf-8')  # 新建该fasta文件，准备存储挑出的潜在cas文件

    for fasta_i in tqdm(selected_Cas_fasta_dict.keys()):  # 前面建立字典的kyes遍历（fasta文件遍历）
        fasta_file_path = os.path.join(predict_result_path, fasta_i)
        # print(fasta_file_path)

        try:
            proteins_obj = Fasta(fasta_file_path)
            os.remove(proteins_obj.faidx.filename+'.fai')  # 删除对应fai文件。
        except:
            continue

        short_names_list = selected_Cas_fasta_dict[fasta_i][0]  # 根据字典的当前键提取对应值的第0个元素（该fasta文件中被选中的蛋白质名称）

        for name_i in short_names_list:
            try:
                name = proteins_obj[name_i].long_name  # 从原始fasta文件中提取对应氨基酸短名称的长名称
                # 写入氨基酸名称
                f.write('>'+name.strip()+'\n') 
                # 写入氨基酸序列
                for seq_seg in proteins_obj[name_i]:  # 逐行提取某个蛋白质的氨基酸序列
                    f.write(str(seq_seg)+'\n')

                num_of_successes +=1
            except:
                print(f"Warning: {name_i} in {fasta_i} 失败！")


    f.close()

print(f'共得到 {num_of_successes} 个潜在Cas1序列！')
