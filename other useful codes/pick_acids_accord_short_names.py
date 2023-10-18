## 说明：从原始fasta文件中，根据‘短名称’选择氨基酸序列，并存储成新的fasta文件。

from pyfaidx import Fasta
import os

# 函数
def pick_acids_accord_seq_names(proteins_obj, new_fasta_file, short_names_list, long_name=True):
    
    '''函数功能：从某个fasta文件中，根据氨基酸序列的“短名称”提取对应氨基酸序列，并存储成新的fasta文件；
        proteins_obj：Fasta(xxx.fasta)读取的fasta原文件对象；
        new_fasta_file：用于存储选中氨基酸序列的文件（存成fasta格式），命名时请以.fasta结尾；
        short_names_list：包含氨基酸序列名称的list, []；
        long_name：新的fasta文件中氨基酸序列信息是按长名称还是短名称存储。
    '''
    
    if not isinstance(short_names_list,list):
        raise ValueError("Warning:'short_names_list' should be a list type.")

    seq_fail = []  # 存放提取失败的氨基酸序列名称

    with open(new_fasta_file, 'w', encoding='utf-8') as f:

        for name_i in short_names_list:
            try:
                if long_name==True:
                    name = proteins_obj[name_i].long_name  # 从原始fasta文件中提取对应氨基酸序列短名称的长名称
                else:
                    name = proteins_obj[name_i].name  # 从原始fasta文件中提取对应氨基酸序列短名称的短名称（实际上就是name_i）
                # print(name)
                
                # 写入氨基酸序列名称
                f.write('>'+name.strip()+'\n')  

                # 写入氨基酸序列
                for seq_seg in proteins_obj[name_i]:  # 逐行提取某个蛋白质的氨基酸序列
                    f.write(str(seq_seg)+'\n')

            except:
                seq_fail.append(name_i)

    print(f'New file has been created! {len(short_names_list)-len(seq_fail)} succeed, {len(seq_fail)} failed!')                  

    return seq_fail



# ## 使用（测试）
if __name__ == '__main__':
    # 氨基酸短名称
    short_names_list = ['WP_011752983.1','AAX71862.1','UniRef50_Q6GZX3','UniRef50_Q197F8']
    os.remove(proteins.filename+'.fai')  # 删除因使用Fasta()方法而自动生成的.fai文件，对项目没啥用。

    proteins = Fasta('test.fasta')
    new_fasta_file = 'select_accord_seq_short_names.FASTA'
    _ = pick_acids_accord_seq_names(proteins, new_fasta_file, short_names_list, long_name=False)