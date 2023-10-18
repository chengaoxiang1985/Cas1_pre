import os
import pandas as pd
from txdpy import get_letter
from tqdm import tqdm
import shutil

# 读取前面保存的存储的所有细菌的蛋白质相关信息文件（csv文件）
file_name = 'df_protein_info.csv'
df = pd.read_csv(file_name, encoding='utf-8')

#  从‘细菌名称’中提取第一个单词，它代表了该细菌属于那个种类
bacteria_kinds = []
for bacteria_name in df['细菌名称']:
#     print(bacteria_name.split(' ')[0:2])
    kind_i = ''.join(get_letter(bacteria_name.split(' ')[0])).lower()  # 提取每个细菌名称的第一个单词（它说明了属于哪种细菌）
    bacteria_kinds.append(kind_i)
#     print(kind_i)

# print(bacteria_kinds)

#  将细菌种类信息插入df中留作后用
df.insert(1,'细菌种类',bacteria_kinds)

#  计算总共有多少种类的细菌
bacteria_sets = list(set(bacteria_kinds))
bacteria_sets.sort()  # 根据字母对单词进行排序（修改原数据）
# print(len(bacteria_sets))
# for i in bacteria_sets:
#     print(i)


# 根据16S rRNA构建的系统发育树，目前发现的古细菌（Archaea ）包括5个界或门
# 泉古菌 （Crenarchaeota）；广古菌(Euryarchaeota)；初生古菌(Korarchaeota)；纳古菌(Nanoarchaeota)；奇古菌(Thaumarchaeota)
# 下述文件中包含了这五大类极其子类古菌（收集得可能不够齐全）

archaeas = []
with open('archaea_list.txt','r',encoding='utf-8') as f:
    for i in f.readlines():
        archaeas.append(i.strip().lower())
# print(archaeas)


archaeas_data_list = []  # 挑选出下载的数据中属于古菌的名称
for i in bacteria_kinds:
    if i in archaeas:
        archaeas_data_list.append(i)
print(len(archaeas_data_list))

df_archaea = df[df['细菌种类'].isin(archaeas_data_list)]  # 根据挑选出的古菌名称，将这些古菌的信息单独抽出
df_archaea


# 生成存储古菌数据的文件
archaeas_file_name = 'archaea_data'
try:
    os.makedirs(archaeas_file_name)  
except:
    pass

# 从已经下载好的数据中，将古细菌挑出并存储至一个文件夹中
for i in tqdm(range(df_archaea.shape[0])):

    data_name_i = df_archaea['蛋白质压缩数据名称'].iloc[i]
    data_path = r"./download_bacteria_protein_data/"+data_name_i[0]+'/'+data_name_i  # 形成数据的完整路径
    # print(data_path)
    
    shutil.copy(data_path, archaeas_file_name)  # 将数据data_path拷贝到archaeas_file_name（前面生成）表示的文件夹中
    
    # try:
    #     shutil.copy(data_path, archaeas_file_name)  # 将数据data_path拷贝到archaeas_file_name（前面生成）表示的文件夹中
    # except:
    #     print(f'拷贝失败：{data_name_i}')