# 利用训练好的模型预测新的蛋白
import os
# os.environ["CUDA_VISIBLE_DEVICES"]="-1"  # 禁止使用GPU，即用cpu跑
os.environ["CUDA_VISIBLE_DEVICES"]="1"  # 使用第0张GPU跑

import sys
sys.path.append(os.path.join(os.getcwd(),'chemprop'))  # 加入chemprop包的路径
sys.path.append(os.path.dirname(os.getcwd()))

import chemprop

print('Program starts......')

# # 1.对model1预测的结果，用model2继续进行预测
# save_model_dir = 'Cas1_with_nocas_smiles_checkpoints_2023-6-26'  # 训练时存储checkpoints的路径
# predict_data_path = 'data_for_pre_cas1_use_model2/potential_cas1_0.9_1_12574_smiles.csv'  # 待预测数据文件路径
# predict_result_path = 'data_for_pre_cas1_use_model2/potential_cas1_0.9_1_12574_smiles_predict.csv'  # 待生成的预测结果路径
# save_npz_path = 'data_for_pre_cas1/potential_cas1_0.9_1_12574_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

# predict_arguments = [
#     '--test_path', predict_data_path,  # 存放待预测数据的路径
#     '--preds_path', predict_result_path,  # 存放预测结果的路径
#     '--checkpoint_dir', save_model_dir
# ]

# predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
# preds = chemprop.train.make_predictions(args=predict_args)


# # 2. 用model2对Cas1和uniref50中抽取的长度在400到1300aa的蛋白进行预测
# save_model_dir = 'Cas1_with_nocas_smiles_checkpoints_2023-6-26'  # 训练时存储checkpoints的路径
# predict_data_path = 'data_for_pre_cas1_use_model2/Cas1_uniref50_len400to1300_smiles.csv'  # 待预测数据文件路径
# predict_result_path = 'data_for_pre_cas1_use_model2/Cas1_uniref50_len400to1300_smiles_predict.csv'  # 待生成的预测结果路径
# save_npz_path = 'data_for_pre_cas1/Cas1_uniref50_len400to1300_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

# predict_arguments = [
#     '--test_path', predict_data_path,  # 存放待预测数据的路径
#     '--preds_path', predict_result_path,  # 存放预测结果的路径
#     '--checkpoint_dir', save_model_dir
# ]

# predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
# preds = chemprop.train.make_predictions(args=predict_args)


# # 3. 用model2对Cas1和部分其他Cas蛋白（Cas2-Cas14）进行预测。
# save_model_dir = 'Cas1_with_nocas_smiles_checkpoints_2023-6-26'  # 训练时存储checkpoints的路径
# predict_data_path = 'data_for_pre_cas1_use_model2/Cas1_and_part_other_Cas_smiles.csv'  # 待预测数据文件路径
# predict_result_path = 'data_for_pre_cas1_use_model2/Cas1_and_part_other_Cas_smiles_predict.csv'  # 待生成的预测结果路径
# save_npz_path = 'data_for_pre_cas1/Cas1_and_part_other_Cas_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

# predict_arguments = [
#     '--test_path', predict_data_path,  # 存放待预测数据的路径
#     '--preds_path', predict_result_path,  # 存放预测结果的路径
#     '--checkpoint_dir', save_model_dir
# ]

# predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
# preds = chemprop.train.make_predictions(args=predict_args)


# 4. 用model2对小于400aa的Cas1和Uniref50进行预测。
save_model_dir = 'Cas1_with_nocas_smiles_checkpoints_2023-6-26'  # 训练时存储checkpoints的路径
predict_data_path = 'data_for_pre_cas1_use_model2/Cas1_Uniref50_len400_smiles.csv'  # 待预测数据文件路径
predict_result_path = 'data_for_pre_cas1_use_model2/Cas1_Uniref50_len400_smiles_predict.csv'  # 待生成的预测结果路径
save_npz_path = 'data_for_pre_cas1/Cas1_Uniref50_len400_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

predict_arguments = [
    '--test_path', predict_data_path,  # 存放待预测数据的路径
    '--preds_path', predict_result_path,  # 存放预测结果的路径
    '--checkpoint_dir', save_model_dir
]

predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
preds = chemprop.train.make_predictions(args=predict_args)