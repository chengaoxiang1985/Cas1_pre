# 跑一个完整的训练模型
import os
# os.environ["CUDA_VISIBLE_DEVICES"]="-1"  # 禁止使用GPU，即用cpu跑
os.environ["CUDA_VISIBLE_DEVICES"]="1"  # 使用第几张GPU（从0开始编号）跑

import sys
sys.path.append(r"/data/home/chengaoxiang/chemprop")

import chemprop

if __name__ == '__main__':

    ## 测试数据
    # arguments = [
    #     '--data_path', 'data/tox21.csv',   # 参数在chemprop-args.py-TrainArg类中
    #     '--dataset_type', 'classification', # 参数在chemprop-args.py-TrainArg类中 dataset_type: Literal['regression', 'classification', 'multiclass', 'spectra']
    #     '--save_dir', 'tox21_checkpoints']
    #
    # # arguments = [
    # #     '--data_path', 'Cas_data/Cas7/Cas7_smiles.csv',  # 参数在chemprop-args.py-TrainArg类中
    # #     '--dataset_type', 'classification',  # 参数在chemprop-args.py-TrainArg类中 dataset_type: Literal['regression', 'classification', 'multiclass', 'spectra']
    # #     '--save_dir', 'Cas7_cgx_checkpoints']
    #
    # args = chemprop.args.TrainArgs().parse_args(arguments)
    # mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)

##################################################################################################################################################
    # # Cas1
    # train_data_path = 'Cas_data/Cas1/Cas1_with_nocas_smiles.csv'  # 训练数据路径
    # dataset_type = 'classification'  # 任务类型
    # log_dir = 'Cas1_with_nocas_smiles_hyperparameter'  # 将写入所有超参数优化结果的目录的路径
    # config_save_path = log_dir + '/' + 'Cas1_with_nocas_smiles_hyperparameter.json'  # 存储超参数的文件，是一个.json文件
    # save_model_dir = 'Cas1_with_nocas_smiles_checkpoints'  # 训练时存储checkpoints的路径

    # # Cas1(positive <400aa) and Uniref50(negative <400aa)
    # train_data_path = 'Cas1_data/Cas1_Uniref50_len400_smiles.csv'  # 训练数据路径
    # dataset_type = 'classification'  # 任务类型
    # log_dir = 'Cas1_with_nocas_smiles_hyperparameter'  # 将写入所有超参数优化结果的目录的路径
    # config_save_path = log_dir + '/' + 'Cas1_with_nocas_smiles_hyperparameter.json'  # 存储超参数的文件，是一个.json文件
    # save_model_dir = 'Cas1_with_nocas_smiles_checkpoints'  # 训练时存储checkpoints的路径

    # Cas1(positive <400aa) and Cas2toCas14(part, negative <600aa)
    train_data_path = 'Cas1_data/Cas1_len400_Cas2toCas14_len600_smiles.csv'  # 训练数据路径
    dataset_type = 'classification'  # 任务类型
    log_dir = 'Cas1_with_nocas_smiles_hyperparameter'  # 将写入所有超参数优化结果的目录的路径
    config_save_path = log_dir + '/' + 'Cas1_with_nocas_smiles_hyperparameter.json'  # 存储超参数的文件，是一个.json文件
    save_model_dir = 'Cas1_with_nocas_smiles_checkpoints'  # 训练时存储checkpoints的路径


    # # Cas8
    # train_data_path = 'Cas_data/Cas8/Cas8_with_nocas_smiles.csv'  # 训练数据路径
    # dataset_type = 'classification'  # 任务类型
    # log_dir = 'Cas8_with_nocas_smiles_hyperparameter'  # 将写入所有超参数优化结果的目录的路径
    # config_save_path = log_dir + '/' + 'Cas8_with_nocas_smiles_hyperparameter.json'  # 存储超参数的文件，是一个.json文件
    # save_model_dir = 'Cas8_with_nocas_smiles_checkpoints'  # 训练时存储checkpoints的路径

    # 设置运行参数
    arguments = [
        '--data_path', train_data_path,  # 参数在chemprop-args.py-TrainArg类中
        '--dataset_type', dataset_type,  # 参数在chemprop-args.py-TrainArg类中 dataset_type: Literal['regression', 'classification', 'multiclass', 'spectra']
        '--save_dir', save_model_dir,
        '--num_folds', '3',
        '--hidden_size', '1200',
        '--depth', '3',
        '--dropout', '0.3',
        '--ensemble_size', '5',
        '--ffn_num_layers', '3',
        '--num_workers', '8',
        '--batch_size', '80',
        '--epochs', '30']

    # 参数解析
    args = chemprop.args.TrainArgs().parse_args(arguments)  # parse_args()是tap自带的函数
    # print(args)

    # 开始训练模型
    mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)
##################################################################################################################################################
