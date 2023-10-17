import os
# os.environ["CUDA_VISIBLE_DEVICES"]="-1"  # 禁止使用GPU，即用cpu跑
os.environ["CUDA_VISIBLE_DEVICES"]="0"  # 使用第0张GPU跑

import sys
sys.path.append(os.path.join(os.getcwd(),'chemprop'))  # 加入chemprop包的路径
sys.path.append(os.path.dirname(os.getcwd()))

import chemprop
import re
from tqdm import tqdm

if __name__=="__main__":
    print('Program starts......')
    
    save_model_dir = 'Cas1_with_nocas_smiles_checkpoints'  # 训练时存储checkpoints的路径
    csv_file_path = '/data/home/chengaoxiang/dload_bacteria_protein_data_cgx/archaea_data/selected_archaea_cgx'  # 待预测的csv文件的地址
    csv_flle_names = [i for i in os.listdir(csv_file_path) if i[-4:]=='.csv']
    for csv_i in tqdm(csv_flle_names):  # 每个csv文件循环
        csv_i_path = os.path.join(csv_file_path,csv_i)  # 待预测文件的完整路径
        # print(csv_i_path)
        result = re.match('^(.*).csv', csv_i)  # 截取文件名部分作为预测结果文件名的一部分
        csv_i_predict_path = os.path.join(csv_file_path, result.group(1)+'.predict.csv')  # 待生成预测结果文件的全路径
        # print(csv_i_predict_path)

        predict_arguments = [
            '--test_path', csv_i_path,  # 存放待预测数据的路径
            '--preds_path', csv_i_predict_path,  # 存放预测结果的路径
            '--checkpoint_dir', save_model_dir  # 训练时存储checkpoints的路径
        ]

        predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
        preds = chemprop.train.make_predictions(args=predict_args)