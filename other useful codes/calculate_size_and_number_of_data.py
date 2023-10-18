# # 计算细菌、古细菌、氨基酸序列个数及大小
import os
import re
from pyfaidx import Fasta
from tqdm import tqdm


if __name__=='__main__':
    print('Program starting......')

    # ## 计算download_bacteria_protein_data路径下所有的蛋白压缩文件的个数及大小（包括了ensemble bacteria中下载的所有细菌和古细菌）
    path = './download_bacteria_protein_data'
    total_file_size = 0
    total_file_num = 0
    for root, dirs, files in tqdm(os.walk(path)):  # 获取路径下的所有子文件夹、子文件
    #     print(root)
    #     print(dirs)
    #     print(files)
        if files:  # 如果获取到子文件
            total_file_num += len(files)
            for file in files:
                # print(file)
                file_size = os.path.getsize(os.path.join(root,file)) 
                # print(file_size)
                total_file_size += file_size

    print(f"共{total_file_num}个细菌和古细菌压缩数据文件，共占{total_file_size/1000000}M字节！")                   



    # ## 计算古细菌的压缩数据文件的个数大小、氨基酸序列总个数
    path = './archaea_data/'

    # 计算压缩古细菌数据的大小和个数
    compressed_archaea_datafile = [file for file in os.listdir(path) if file[-3:]=='.gz']  # 提取压缩文件名称
    toal_comp_num = len(compressed_archaea_datafile)  # 压缩文件个数
    total_comp_size = 0  # 总大小
    for file in tqdm(compressed_archaea_datafile):
        file_size = os.path.getsize(os.path.join(path,file))
        total_comp_size +=file_size
    print(f"共{toal_comp_num}个压缩古细菌数据，共占{total_comp_size/1000000}M字节！")

    # 计算解压的古细菌数据的大小和个数，并统计所有氨基酸序列个数
    uncompressed_archaea_datafile = [file for file in os.listdir(path) if file[-3:]=='.fa']  # 提取非压缩文件名称
    toal_uncomp_num = len(uncompressed_archaea_datafile)  # 非压缩文件个数
    total_uncomp_size = 0  # 总大小
    total_protein_num = 0  # 总氨基酸序列个数
    for file in tqdm(uncompressed_archaea_datafile):
        file_size = os.path.getsize(os.path.join(path,file))
        total_uncomp_size +=file_size
        protein_obj = Fasta(os.path.join(path,file))  # 打开对应的fasta文件
        total_protein_num +=len(protein_obj.keys())  # 计算该文件中共有多少个氨基酸序列
        os.remove(protein_obj.filename+'.fai')  # 删除因调用Fasta()而产生的.fai文件
        
    print(f"共{toal_uncomp_num}个解压古细菌数据，共占{total_uncomp_size/1000000}M字节！")
    print(f"共{total_protein_num}个古细菌氨基酸序列！")



    # ## 计算被选出的古细菌的氨基酸序列个数（用于寻找Cas1）
    path = './archaea_data/selected_archaea_cgx/'
    selected_archaea_fasta = [file for file in os.listdir(path) if file[-16:]=='.cgxselect.fasta']  # 提取压缩文件名称
    total_protein_num = 0  # 总氨基酸序列个数
    for file in tqdm(selected_archaea_fasta):
        file_size = os.path.getsize(os.path.join(path,file))
        protein_obj = Fasta(os.path.join(path,file))  # 打开对应的fasta文件
        total_protein_num +=len(protein_obj.keys())  # 计算该文件中共有多少个氨基酸序列
        os.remove(protein_obj.filename+'.fai')  # 删除因调用Fasta()而产生的.fai文件

    print(f"共{total_protein_num}个被选中的古细菌氨基酸序列（用于搜索是否为Cas1）！")