# Cas1_pre
A graph neural network model for predicting Cas1 protein.

Description
-----------
![image](https://github.com/chengaoxiang1985/Cas1_pre/blob/main/graph%20abstract.png)
The CRISPR-Cas system, an adaptive immune mechanism found in bacteria and archaea, has evolved into a promising genome editing tool. Cas proteins, including Cas1, play vital roles in acquiring spacer sequences and integrating foreign nucleic acids. In this study, we first gathered and analyzed a comprehensive collection of CRISPR-associated (Cas) proteins, ranging from Cas1 to Cas14. Specifically, we focused on Cas1 and converted these proteins into the simplified molecular-input line-entry system (SMILES) format to construct graph data representing atom and bond features. Next, two GNN models were designed using the directed message passing neural network (DMPNN) framework, and these models were trained on two carefully curated Cas1 graph datasets. Subsequently, the performance of these models on both the training data and newly designed datasets was evaluated, and then compared with a widely used non-deep learning method. Finally, the established models were used to identify new Cas1 proteins within the Ensemble database. Our models demonstrated their effectiveness in identifying previously unknown Cas1 proteins, highlighting their robustness and practical utility. In conclusion, our models serve as a valuable auxiliary tool for Cas1 protein identification, and contribute to the innovative application of SMILES encoding in the study of biomacromolecules.

Raw and training data
----------
**Cas_data** contains all the Cas proteins we have collected, from Cas1 to Cas14; **Cas_data_TSNE** performs tsne analysis based on **Cas_data**; **Cas1_data** contains all the training data used for our project.  **Bacteria protein URL** is where to download bacteria proteins.  
- **Cas_data:** https://drive.google.com/drive/folders/1ZU3BeMpCQ15VLfPMczuYdwQVjSKHlIAB?usp=drive_link  
- **Cas_data_TSNE:** https://drive.google.com/drive/folders/1hr-kXUi5O02oUz5gAFzL2p5UBT2sx_zW?usp=drive_link  
- **Cas1_data:** https://drive.google.com/drive/folders/1beRpDvDpfurzGWFecagMysnScKnh0Vw-?usp=drive_link  
- **Bacteria protein URL:** https://drive.google.com/file/d/1lYm5M0ODwrpyVjmg5m0Lg8UxolwzYDYs/view?usp=share_link

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
1.Install all packages listed on this page.  
2.Download or clone the repository.
3.Follow the examples on this page to retain new model or use our trained models.
(It is important to note that some paths in the code file need to be adjusted according to the actual situation.)

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
