# Cas1_pre
A graph neural network model for predicting Cas1 protein.


Description
-----------
![image](https://github.com/chengaoxiang1985/Cas1_pre/blob/main/graph%20abstract.png)

Raw and training data
----------
**Cas_data** contains all the Cas proteins we have collected, from Cas1 to Cas14; **Cas_data_TSNE** performs tsne analysis based on **Cas_data**; **Cas1_data** contains all the training data used for our project.  
- **Cas_data:** https://drive.google.com/drive/folders/1ZU3BeMpCQ15VLfPMczuYdwQVjSKHlIAB?usp=drive_link  
- **Cas_data_TSNE:** https://drive.google.com/drive/folders/1hr-kXUi5O02oUz5gAFzL2p5UBT2sx_zW?usp=drive_link  
- **Cas1_data:** https://drive.google.com/drive/folders/1beRpDvDpfurzGWFecagMysnScKnh0Vw-?usp=drive_link  

Data for prediction
----------
- **data_for_pre_cas1** the data set predicted by model 1 and the predicted results; **data_for_pre_cas1_use_model2** the data set predicted by model 2 and the predicted results.  
- **data_for_pre_cas1:** https://drive.google.com/drive/folders/18oOQ4n5EGcr8EYFNUOvKheCLtSuYb0aq?usp=drive_link  
- **data_for_pre_cas1_use_model2:** https://drive.google.com/drive/folders/16BSfr1AxMgMd3HXE3f6mQXeUGqCF3W21?usp=drive_link  

Moldel checkpoints
----------
- **Model 1: Cas1_with_nocas_smiles_checkpoints_2023-5-8**  
  https://drive.google.com/drive/folders/1UFxpk8QgewgIrr159hCLDCb3eSFpFQx0?usp=drive_link  
- **Model 2: Cas1_with_nocas_smiles_checkpoints_2023-6-26**  
  https://drive.google.com/drive/folders/1mF9oA4wefPXFKAYFhnVi-V5f4cTrN2tJ?usp=drive_link  

Training example
----------
```python
import os
# os.environ["CUDA_VISIBLE_DEVICES"]="-1"  
os.environ["CUDA_VISIBLE_DEVICES"]="1"  

import sys
sys.path.append(r"/data/home/chengaoxiang/chemprop")

import chemprop

if __name__ == '__main__':
    arguments = [
        '--data_path', train_data_path, 
        '--dataset_type', dataset_type, 
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

    args = chemprop.args.TrainArgs().parse_args(arguments) 
    # print(args)

    # training start
    mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)
```

Predicting example
----------
```python

import os
# os.environ["CUDA_VISIBLE_DEVICES"]="-1"
os.environ["CUDA_VISIBLE_DEVICES"]="1"

import sys
sys.path.append(os.path.join(os.getcwd(),'chemprop')) 
sys.path.append(os.path.dirname(os.getcwd()))

import chemprop

print('Program starts......')
predict_arguments = [
    '--test_path', predict_data_path,
    '--preds_path', predict_result_path,
    '--checkpoint_dir', save_model_dir
]

predict_args = chemprop.args.PredictArgs().parse_args(predict_arguments)
preds = chemprop.train.make_predictions(args=predict_args)
```

Usage
----------
** 1.fsw  
2.sdfa  
3.asdf  



Main packages used
----------
- Python 3.9.13  
- Pytorch 1.12.1.  
- chemprop
- rdkit
- numpy
- csv
- tqdm

The hardware and os we used
----------
- Intel® Xeon(R) Gold 6248R CPU @ 3.00GHz × 4 (Intel, Santa Clara, CA, USA)  
- Ubuntu 7.5.0-3ubuntu1~18.04 (Linux version 4.15.0-112-generic) and 754 GB RAM  
- Four NVIDIA Tesla V100S PCIe GPUs (32GB each, CUDA Version: 11.6)

Acknowledgements
----------
This research is supported by Key Research Project of Zhejiang Lab (No. 117005-AC2106/001)
