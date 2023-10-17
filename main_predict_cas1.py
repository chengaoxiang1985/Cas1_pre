# 利用训练好的模型预测新的蛋白
import os
# os.environ["CUDA_VISIBLE_DEVICES"]="-1"  # 禁止使用GPU，即用cpu跑
os.environ["CUDA_VISIBLE_DEVICES"]="0"  # 使用第0张GPU跑

import sys
sys.path.append(os.path.join(os.getcwd(),'chemprop'))  # 加入chemprop包的路径
sys.path.append(os.path.dirname(os.getcwd()))

import chemprop

print('Program starts......')

# # 1.预测真假Cas1，即训练数据（长度都是400）
# save_model_dir = 'Cas1_with_nocas_smiles_checkpoints'  # 训练时存储checkpoints的路径
# predict_data_path = 'data_for_pre_cas1/Cas1_Uniref50_len400_smiles.csv'  # 待预测数据文件路径
# predict_result_path = 'data_for_pre_cas1/Cas1_Uniref50_len400_smiles_predict.csv'  # 待生成的预测结果路径
# save_npz_path = 'data_for_pre_cas1/Cas1_Uniref50_len400_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

# predict_arguments = [
#     '--test_path', predict_data_path,  # 存放待预测数据的路径
#     '--preds_path', predict_result_path,  # 存放预测结果的路径
#     '--checkpoint_dir', save_model_dir
# ]

# predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
# preds = chemprop.train.make_predictions(args=predict_args)


# # 2.预测真假Cas1，即训练数据（氨基酸长度1-800）
# save_model_dir = 'Cas1_with_nocas_smiles_checkpoints'  # 训练时存储checkpoints的路径
# predict_data_path = 'data_for_pre_cas1/Cas1_Uniref50_len800_smiles.csv'  # 待预测数据文件路径
# predict_result_path = 'data_for_pre_cas1/Cas1_Uniref50_len800_smiles_predict.csv'  # 待生成的预测结果路径
# save_npz_path = 'data_for_pre_cas1/Cas1_Uniref50_len800_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

# predict_arguments = [
#     '--test_path', predict_data_path,  # 存放待预测数据的路径
#     '--preds_path', predict_result_path,  # 存放预测结果的路径
#     '--checkpoint_dir', save_model_dir
# ]

# predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
# preds = chemprop.train.make_predictions(args=predict_args)

# # 3.预测真假Cas1，非训练数据（氨基酸长度800以上）
# save_model_dir = 'Cas1_with_nocas_smiles_checkpoints'  # 训练时存储checkpoints的路径
# predict_data_path = 'data_for_pre_cas1/Cas1_uniref50_len400to1300_smiles.csv'  # 待预测数据文件路径
# predict_result_path = 'data_for_pre_cas1/Cas1_uniref50_len400to1300_smiles_predict.csv'  # 待生成的预测结果路径
# save_npz_path = 'data_for_pre_cas1/Cas1_uniref50_len400to1300_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

# predict_arguments = [
#     '--test_path', predict_data_path,  # 存放待预测数据的路径
#     '--preds_path', predict_result_path,  # 存放预测结果的路径
#     '--checkpoint_dir', save_model_dir
# ]

# predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
# preds = chemprop.train.make_predictions(args=predict_args)

# # 4.预测真假Cas1，(真Cas1是所有长度的Cas1，假Cas1是从所有其他每个非Cas1（Cas2-Cas14）类型中抽出的部分序列)
# save_model_dir = 'Cas1_with_nocas_smiles_checkpoints'  # 训练时存储checkpoints的路径
# predict_data_path = 'data_for_pre_cas1/Cas1_and_part_other_Cas_smiles.csv'  # 待预测数据文件路径
# predict_result_path = 'data_for_pre_cas1/Cas1_and_part_other_Cas_smiles_predict.csv'  # 待生成的预测结果路径
# save_npz_path = 'data_for_pre_cas1/Cas1_and_part_other_Cas_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

# predict_arguments = [
#     '--test_path', predict_data_path,  # 存放待预测数据的路径
#     '--preds_path', predict_result_path,  # 存放预测结果的路径
#     '--checkpoint_dir', save_model_dir
# ]

# predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
# preds = chemprop.train.make_predictions(args=predict_args)


# # 5.预测真Cas1，(真Cas1是被选择用于训练模型的Cas1，非全部Cas1)
# save_model_dir = 'Cas1_with_nocas_smiles_checkpoints_2023-5-8'  # 训练时存储checkpoints的路径
# predict_data_path = 'data_for_pre_cas1/Cas1_len400_train_smiles.csv'  # 待预测数据文件路径
# predict_result_path = 'data_for_pre_cas1/Cas1_len400_train_smiles_predict.csv'  # 待生成的预测结果路径
# save_npz_path = 'data_for_pre_cas1/Cas1_len400_train_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

# predict_arguments = [
#     '--test_path', predict_data_path,  # 存放待预测数据的路径
#     '--preds_path', predict_result_path,  # 存放预测结果的路径
#     '--checkpoint_dir', save_model_dir
# ]

# predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
# preds = chemprop.train.make_predictions(args=predict_args)



# # 6.预测真Cas1，(真Cas1，长度大于800aa)
# save_model_dir = 'Cas1_with_nocas_smiles_checkpoints_2023-5-8'  # 训练时存储checkpoints的路径
# predict_data_path = 'data_for_pre_cas1/Cas1_len800_smiles.csv'  # 待预测数据文件路径
# predict_result_path = 'data_for_pre_cas1/Cas1_len800_smiles_predict.csv'  # 待生成的预测结果路径
# save_npz_path = 'data_for_pre_cas1/Cas1_len800_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

# predict_arguments = [
#     '--test_path', predict_data_path,  # 存放待预测数据的路径
#     '--preds_path', predict_result_path,  # 存放预测结果的路径
#     '--checkpoint_dir', save_model_dir
# ]

# predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
# preds = chemprop.train.make_predictions(args=predict_args)


# # 7.预测真Cas1，(真Cas1，长度大于800aa)
# save_model_dir = 'Cas1_with_nocas_smiles_checkpoints_2023-5-8'  # 训练时存储checkpoints的路径
# predict_data_path = 'data_for_pre_cas1/Cas1_uniref50_len400to1300_smiles.csv'  # 待预测数据文件路径
# predict_result_path = 'data_for_pre_cas1/Cas1_uniref50_len400to1300_smiles_predict.csv'  # 待生成的预测结果路径
# save_npz_path = 'data_for_pre_cas1/Cas1_uniref50_len400to1300_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

# predict_arguments = [
#     '--test_path', predict_data_path,  # 存放待预测数据的路径
#     '--preds_path', predict_result_path,  # 存放预测结果的路径
#     '--checkpoint_dir', save_model_dir
# ]

# predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
# preds = chemprop.train.make_predictions(args=predict_args)


# 8.预测真Cas1，(真Cas1，长度大于800aa)
save_model_dir = 'Cas1_with_nocas_smiles_checkpoints_2023-5-8'  # 训练时存储checkpoints的路径
predict_data_path = 'data_for_pre_cas1/potential_cas1_0.9_1_12574_smiles.csv'  # 待预测数据文件路径
predict_result_path = 'data_for_pre_cas1/potential_cas1_0.9_1_12574_smiles_predict.csv'  # 待生成的预测结果路径
save_npz_path = 'data_for_pre_cas1/potential_cas1_0.9_1_12574_smiles_features.npz'   # Path to .npz file where features will be saved as a compressed numpy archive(待生成文件名和路径)

predict_arguments = [
    '--test_path', predict_data_path,  # 存放待预测数据的路径
    '--preds_path', predict_result_path,  # 存放预测结果的路径
    '--checkpoint_dir', save_model_dir
]

predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
preds = chemprop.train.make_predictions(args=predict_args)