#!/usr/bin/env python
# coding: utf-8

# ## 说明：将一个较大的fasta文件，按顺序均匀分割成多个较小的fasta文件

# ## 定义函数

import time
import numpy as np
import os


def split_fasta_larg2small(old_fasta_name:str, new_fasta_name_prefix:str, new_fasta_name_suffix:str, size:int, new_data_file='new_data_file'):
    '''
    该函数将一个较大的fasta文件，按顺序均匀分割成多个较小的fasta文件
    old_fasta_name：原fasta文件名
    new_fasta_name_prefix：新fasta文件名前缀
    new_fasta_name_suffix：新fasta文件名后缀（通常设置为：.fasta）(前缀和后缀之间将以阿拉伯数字作为区分)
    size：单个新fasta文件中计划存储多少个氨基酸序列（如果该值比原fasta文件中所有的氨基酸总和还要大，则把所有氨基酸写入一个新的fasta文件中）
    new_data_file：make a new file for new data in current file.
    '''

    # make a new file for new data in current file.
    os.mkdir(new_data_file)

    # 判断输入size是否满足要求
    if not (isinstance(size,int) and size>0):
        raise Exception('Input error: \'size\' must be int and large than 0!')

    file_count = 0  # 记录总共生成了多少个fasta文件
    total_seq_count = 0  # 记录原fasta文件中总共有多少个氨基酸序列
    
    with open(old_fasta_name, encoding='utf-8') as old_f:
        seq_count = 0  # 记录一个新fasta文件中已经写入了多少个氨基酸
        for line in old_f:
            # time.sleep(0.000000000000000000001)
            if len(line.strip()) > 0:  # 确保当前处理行不是空行
                if line.strip()[0] == ">":  # 如果当前行是说明行（说明行的末尾带有“>”符号）
                    if (int(np.floor(seq_count/size)) == 1 or seq_count == 0):  # 判断是否需要换一个文件开始存储

                        try: # 如果存在上一个数据文件，则应先关闭
                            new_f.close()
                        except:
                            pass

                        new_file_name = new_data_file + '/' + new_fasta_name_prefix+'_'+str(file_count)+new_fasta_name_suffix  # 以数字作为分割的数据文件名
                        new_f = open(new_file_name, 'w', encoding='utf-8')
                        new_f.write(line.strip())  # 将说明行写入文件
                        new_f.write('\n')  # 说明行换行（否则和氨基酸序列连一起）
                        file_count += 1  # 新文件计算变量+1
                        seq_count = 0  # 新产生一个fasta文件时，新文件内氨基酸计数变量清零

                    else:  # 如果不需要换文件存储，则继续将该行存储在打开的文件中
                        try:
                            new_f.write(line.strip())  # 将说明行写入文件
                            new_f.write('\n')  # 说明行换行（否则和氨基酸序列连一起）
                        except:
                            pass

                    seq_count += 1  # 新文件内氨基酸计数变量+1
                    total_seq_count += 1  # 总氨基酸计数变量+1

                else:  # 当前行是氨基酸序列
                    try:
                        new_f.write(line.strip())  # 将氨基酸信息写入新文件
                        new_f.write('\n') # 氨基酸换行存储（标准的fasta文件要求单行的序列长度不超过80个字符）
                                          #（该句如果注释掉，则每个氨基酸序列将存成1行）★★★★★
                    except:
                        pass
            else:
                print('注意：有空值！')

        new_f.close()  # 关闭最后一个数据文件
        
    return file_count, total_seq_count


def run_time(seconds:int):
    '''
    将seconds（秒数）转换为秒、分、小时，并输出信息
    seconds：秒数
    '''
    # 判断输入seconds是否满足要求
    if not (seconds>0):
        raise Exception('Input error: \'seconds\' must large than or equal to 0!')
    
    if seconds<60:
        print(f'共耗时{seconds}秒！')
    if 60<=seconds<3600:
        print(f'共耗时{seconds/60}分钟！')
    if 3600<=seconds:
        print(f'共耗时{seconds/3600}小时！')


# ## 使用

if __name__ == '__main__':
    print('start...')

    # # 小型fasta文件测试
    # start_t = time.time()  # 起始时间
    #
    # file_count, total_seq_count = split_fasta_larg2small(old_fasta_name='test2.fasta',  # 源文件
    #                                                      new_fasta_name_prefix='cgx',  # 新文件前缀
    #                                                      new_fasta_name_suffix='.fasta',  # 新文件后缀
    #                                                      new_data_file='new_data_file',  # 新文件后缀
    #                                                      size=3)  # 单个新文件包含的氨基酸序列个数
    # print(f'共生成{file_count}个新fasta文件！')
    # print(f'原fasta文件共包含{total_seq_count}个氨基酸！')
    #
    # end_t = time.time()  # 终止时间
    # run_time(end_t-start_t)

    # 大fasta文件测试
    # 对uniref50.fasta数据文件（大小为20G左右，包含大约54150000个氨基酸数据）进行切分，安顺序将10000个氨基酸序列存于一个fasta文件中，共生成5415个独立文件。
    # 注意因为文件较大，执行时间比较长
    start_t = time.time()  # 起始时间

    file_count, total_seq_count = split_fasta_larg2small(old_fasta_name='Cas1_combine.FASTA',  # 源文件
                                                         new_fasta_name_prefix='uniref50',  # 新文件前缀
                                                         new_fasta_name_suffix='.fasta',  # 新文件后缀
                                                         new_data_file='uniref50_len800_splited',  # 新建文件夹，存储数据
                                                         size=100)  # 单个新文件包含的氨基酸序列个数
    print(f'共生成{file_count}个新fasta文件！')
    print(f'原fasta文件共包含{total_seq_count}个氨基酸！')

    end_t = time.time()  # 终止时间
    run_time(end_t-start_t)