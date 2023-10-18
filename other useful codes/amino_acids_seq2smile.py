import os
from rdkit import Chem
import csv
from tqdm import tqdm
from pyfaidx import Fasta


# ## 将氨基酸序列数据转化成smiles数据¶
# - 程序将每个氨基酸序列转换成了对应的smiles格式数据

# 数据文件路径
path = os.path.join(os.getcwd(),'uniref50_len800_splited')
print(path)

most_seqlen_2_smiles = 500  #★★★★★★★★★(经测试，在本机中，氨基酸序列大于600时，用rdkit转化为smiles将导致内核奔溃)

seq_len_tongji = []  # 存放每个氨基酸序列的长度
smiles_len_tongji = []  # 存放每个smiles序列的长度
num_seq2smi = 0  # 成功被转化成smiles序列的氨基酸

for file_i in tqdm(os.listdir(path)):
    #     print(file_i, flush=True)
    if file_i[-6:] == '.fasta':  # 判断是否为fasta文件（文件夹中还包括其他文件）
        seq_smiles = os.path.join(path,f'{file_i[:-6]}_smiles.csv')  # 为某fasta文件生成对应的csv文件

        with open(seq_smiles, "w", newline="") as datacsv:   # 新建csv文件
            # dialect为打开csv文件的方式，默认是excel，delimiter="\t"参数指写入的时候的分隔符
            csvwriter = csv.writer(datacsv, dialect=("excel"))
            # #写入第一行信息，即列名称。
            csvwriter.writerow(["smiles", "amino_acides_name"])

            F = Fasta(os.path.join(path,file_i))
            for key in F.keys():  # fasta文件中所有氨基酸名称循环
                seq = str(F[key][:])  # 转化成字符串，否则Chem.MolFromFASTA(seq)无法识别
                #print(seq)   # 根据氨基酸名称挨个读取氨基酸序列
                
                seq_len_tongji.append(len(seq)) 

                if 0 < len(seq) < most_seqlen_2_smiles:  # ★★如果序列大于600将出现内核奔溃★★★★★★★★★★★★★★
                    mol = Chem.MolFromFASTA(seq)  # 将序列转化成smiles格式数据
                    if mol is not None:
                        seq_smiles = Chem.MolToSmiles(mol)
                        num_seq2smi += 1  # 记录转换成功的
                        smiles_len_tongji.append(len(seq_smiles))  # 记录当前smiles序列长度
                        # 将smiles格式数据保存到前面创建的csv文件中,并同时指定标签1
                        csvwriter.writerow([seq_smiles, key])
                        
print(f'总共{len(seq_len_tongji)}个氨基酸序列，转换成功{num_seq2smi}个！')
print(f'氨基酸最长长度：{max(seq_len_tongji)}')
print(f'smiles最长长度：{max(smiles_len_tongji)}')


from matplotlib import pyplot as plt

# 绘图
fig,ax = plt.subplots(1,1)
plt.hist(seq_len_tongji, bins=20)
ax.set_xlabel('seq index', fontdict=dict(fontsize=16, color='r'))
ax.set_ylabel('seq length', fontdict=dict(fontsize=16, color='r'))
plt.tick_params(axis='x', colors='red' )
plt.tick_params(axis='y', colors='red' )
plt.show()

plt.hist(smiles_len_tongji, bins=20)
plt.xlabel('smiles length', fontdict=dict(fontsize=16, color='r'))
plt.ylabel('number of smiles', fontdict=dict(fontsize=16, color='r'))
plt.tick_params(axis='x' ,colors= 'red' )
plt.tick_params(axis='y' ,colors= 'red' )

plt.show()