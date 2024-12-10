#!/usr/bin/env python
# coding: utf-8

# # 说明：氨基酸序列embedding、然后tSNE降维、然后聚类分析
# - 1）利用Facebook预训练模型'facebook/esm2_t33_650M_UR50D' 对fasta文件内的氨基酸序列进行训练，将氨基酸embedding成数据序列
# - 2）将embedding后的高维数据通过tSNE方法降维成2维数据
# - 3）对降维后的数据进行聚类分析（cluster）
# 
# 注意：esm2_t33_650M_UR50D是facebook专门针对蛋白质序列训练的模型（https://huggingface.co/facebook，该网址可下载Facebook预训练模型）

# # 利用ESM预训练模型embedding氨基酸序列
# - 包括训练，并将结果保存硬盘以便作为他用

# ## 定义函数
import torch
from pyfaidx import Fasta
from transformers import EsmTokenizer, EsmModel
import pandas as pd
import pickle
from tqdm import tqdm
import numpy as np
from sklearn.manifold import TSNE
import plotly.express as px
import plotly.io as pio
from nltk.cluster import KMeansClusterer, util


def build_esm_model(model_name):
    """自动从网站上下载Facebook的预训练EMS模型
    model_name：实际上就是‘facebook/esm2_t33_650M_UR50D’（该预训练模型在网站上的固定名称）
    注意：第一次执行该函数时，会下载相关数据和文件，耗时较大
    耗时主要用于下载huggingface文件夹（3G左右），，
    pc上存储路径：C:\\Users\chengaoxiang\.cache\huggingface\...
    linux server上：/data/home/chengaoxiang/.cache/huggingface/...
    ★★(下载过程可能被中断，中断后需要手动把这整个huggingface文件删除，否则重新执行代码会出错！！）
    """

    tokenizer = EsmTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name)
    return tokenizer, model

def generate_embeddings(tokenizer, model, seqs, max_seq_len, batch_size=10, set_device='cpu'):
    """在预训练模型上继续训练
    seqs:按行存储的氨基酸序列，一行表示一个完整的氨基酸
    set_device: 'cuda' 或 'cpu'
    all_embeddings:训练结束后，与seqs对应的每个氨基酸的embedding值
    """

    device = torch.device(set_device)
    model = model.to(device)
    ids = tokenizer.batch_encode_plus([str(seq) for seq in seqs],
            add_special_tokens=True, padding="longest")
    input_ids = torch.tensor(ids['input_ids']).to(device)
    attention_mask = torch.tensor(ids['attention_mask']).to(device)

    all_embeddings = np.zeros((0, max_seq_len))
    with torch.no_grad():
        for i in tqdm(range(0, len(input_ids), batch_size)):
            batch_ids = input_ids[i:i+batch_size]
            batch_mask = attention_mask[i:i+batch_size]
            outputs = model(input_ids=batch_ids, attention_mask=batch_mask)
            embeddings = outputs.last_hidden_state[:, 0, :].detach().cpu().numpy()
            all_embeddings = np.concatenate((all_embeddings, embeddings), axis=0)
    return all_embeddings

# 聚类函数
def make_clusters(samples, cluster_count):
    '''对数据进行聚类
        samples:样本集，一行一个样本
        cluster_count: 用户指定聚多少类
        assigned_clusters：list格式返回每个样本对应的类型(类型以数字作为标记：0，1，2……)'''

    kclusterer = KMeansClusterer(
        cluster_count, distance=util.euclidean_distance,
        repeats=25, avoid_empty_clusters=True)
    assigned_clusters = kclusterer.cluster(samples, assign_clusters=True)
    return assigned_clusters

def extract_seqs(fasta_name):
    """从fasta文件中抽取出所有氨基酸的名称和序列，并存储在dict中
    cateorys：设计用于标记每个氨基酸的类别，该参数可有可无，也可后续增加
    """
    
    seqs_dict = {'names': [], 'seqs': [], 'categorys': []}
    for seq in Fasta(fasta_name):
        seqs_dict['names'].append(seq.name)
        seqs_dict['seqs'].append(str(seq))
        seqs_dict['categorys'].append('cgx')  # 随便给个类，实际中根据情况而定
    return seqs_dict


# ## 执行代码
if __name__== "__main__" :
    print('Program start, please wait......')

    # # 数据准备########################################################################################################
    # # choice 1 测试数据
    # seqs_dict1 = extract_seqs('cgx1.fasta')  # 读取数据
    # seqs_dict1['categorys'] = len(seqs_dict1['categorys']) * ['Cas1']  # 将类别修改成‘Cas1’
    # seqs_dict2 = extract_seqs('cgx2.fasta')  # 读取数据
    # seqs_dict2['categorys'] = len(seqs_dict2['categorys']) * ['non Cas1']  # 将类别修改成‘non Cas1’

    # # choice 2 实际数据(氨基酸序列长度最高800)
    # seqs_dict1 = extract_seqs('Cas1_len800.FASTA')  # 读取数据
    # seqs_dict1['categorys'] = len(seqs_dict1['categorys'])*['Cas1']  # 将类别修改成‘Cas1’
    # seqs_dict2 = extract_seqs('picked_f_uniref50_len800_uniform.FASTA') # 读取数据
    # seqs_dict2['categorys'] = len(seqs_dict2['categorys'])*['non Cas1']  # 将类别修改成‘non Cas1’

    # choice 3 实际数据(氨基酸序列长度最高400)
    seqs_dict1 = extract_seqs('Cas1_len400_train.FASTA')  # 读取数据
    seqs_dict1['categorys'] = len(seqs_dict1['categorys'])*['Cas1']  # 将类别修改成‘Cas1’
    seqs_dict2 = extract_seqs('picked_f_uniref50_len400_uniform_train.FASTA') # 读取数据
    seqs_dict2['categorys'] = len(seqs_dict2['categorys'])*['non Cas1']  # 将类别修改成‘non Cas1’

    seqs_dict = {}  # 把上述两个数据融合
    seqs_dict['names'] = seqs_dict1['names'] + seqs_dict2['names']
    seqs_dict['seqs'] = seqs_dict1['seqs'] + seqs_dict2['seqs']
    seqs_dict['categorys'] = seqs_dict1['categorys'] + seqs_dict2['categorys']

    seqs = seqs_dict['seqs']  # 提取所有氨基酸序列
    print(len(seqs))


    # # 开始EMS模型训练##################################################################################################
    # 导入预训练模型（第一次执行需要从网上下载，耗时较长）
    model_name = 'facebook/esm2_t33_650M_UR50D'  # facebook专门针对蛋白质序列训练的模型：https://huggingface.co/facebook
    tokenizer, model = build_esm_model(model_name)

    # 开始训练
    max_seq_len = 1280  # 对于'facebook/esm2_t33_650M_UR50D'模型，长度必须是1280
    embeddings = generate_embeddings(tokenizer, model, seqs, max_seq_len, batch_size=25, set_device='cuda')  # for linux server
    # embeddings = generate_embeddings(tokenizer, model, seqs, max_seq_len, batch_size=5, set_device='cpu')  # for windows

    # 保存数据
    seqs_dict['embeddings'] = embeddings.tolist()  # 将embedding后的数据加入原来的字典
    df = pd.DataFrame(seqs_dict)  # 将seqs_dict转化成Dataframe格式，便于保存和绘图
    pickle.dump(df, open('seqs_embdding_from_EMSmodel.bin', 'wb'))  # 保存数据到硬盘以被后用★★
    df.head()


    # # tSNE降维########################################################################################################
    # - 这部分代码和前面无关
    # 导入数据(前面保存在磁盘上)
    df = pickle.load(open('seqs_embdding_from_EMSmodel.bin', 'rb'))
    df.head()
    embeddings = np.array(df['embeddings'].tolist())  # 抽出每个氨基酸序列的embedding值并转化成二维矩阵
    tsne = TSNE(n_components=2, verbose=1, perplexity=30, n_iter=1000, random_state=2)
    tsne_results = tsne.fit_transform(embeddings)
    df['x'] = tsne_results[:, 0] # 将tsne后的结果也存储在data数据中
    df['y'] = tsne_results[:, 1]
    df.head()
    # 此处也可将带有tsne结果的数据df保存至磁盘留做他用

    # 绘图（注意，利用px绘图，关闭notebook后图片消失，必须重新执行程序绘图）
    fig = px.scatter(df, x='x', y='y', color='categorys', symbol="categorys", title='tSNE降维结果')
    fig.show()
    pio.write_html(fig, "tSNE_result.html")  # 将整个图保存为html，打开有仍然是互动图(plotly express 保存成jpg格式的图片比较复杂)。

    # # # kmean聚类########################################################################################################
    # # - 这部分代码和前面无关（但前面没有把tsne后的结果保存到磁盘，因此数据要用到前的）
    # cluster_count = 3  # 指定几类
    # assigned_clusters = make_clusters(tsne_results, cluster_count=cluster_count)  # 聚类
    # df['clusters'] = pd.Series([str(i) for i in assigned_clusters], index=df.index)  # 将聚类结果给出的每个序列的类别也写入df中保存
    # df.head()

    # # 绘图（注意，利用px绘图，关闭notebook后图片消失，必须重新执行程序绘图）
    # fig1 = px.scatter(df, x='x', y='y', color='clusters', symbol="clusters", title='对tSNE降维结果进行聚类的结果')
    # fig1.show()
    # pio.write_html(fig1, "cluster_on_tSNE_result.html")  # 将整个图保存为html，打开有仍然是互动图(plotly express 保存成jpg格式的图片比较复杂)。