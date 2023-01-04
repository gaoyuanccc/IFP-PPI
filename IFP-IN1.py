####
# 时间：20220220
# 地点：华中农业大学
# 姓名：高远
####

####导入系统包
import os
import time
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy import stats
import argparse
#计算时间
start = time.time()


####导入自定义包
#从package文件夹中导入
from package import cid_divide_package
from package import genbank_package
from package import gene_distance_package
from package import cid_total_package
from package import bin_interaction_transiform_gene_package


####参数设置
parser = argparse.ArgumentParser('Gene association network obtained by three-dimensional spatial interaction frequency of bacterial chromosomes.')
parser.add_argument('-i','--input', type=str, default= 'input',help='Input folder name.')
parser.add_argument('-gb','--genbank', type=str, default='' ,help='Genbank filename with annotation information for the strain.')
parser.add_argument('-im','--IFmatrix', type=str,required=True ,help='Contact matrix filename to convert interaction frequencies between bins to interaction frequencies between genes.')
parser.add_argument('-ir','--IFmatrixResolution', type=int,required=True,help='Resolution of IFmatrix.')
parser.add_argument('-b','--binBeginNumber', type=int, default=0 ,help='Begin the bin serial number in the contact matrix, which is 0 by default.')
parser.add_argument('-n','--genepairNumber', type=int, default=10000 ,help='The number of selected gene pairs and is 10000 by default.')
parser.add_argument('-d','--removeDistance', type=int, default=300 ,help='Linear genomic distances to select removed gene pairs and is 300 by default.')
parser.add_argument('-cm','--CIDmatrix', type=str, default= '',help='Contact matrix filename to identify CIDs and is IFmatrix by default. If CIDmatrix = None, then it will not identify CIDs')
parser.add_argument('-cr','--CIDmatrixResolution',type=int,help='Resolution of CIDmatrix and is IFmatrixResolution by default.')
parser.add_argument('-o','--output', type=str, default= 'output',help='Output folder name.')
parser.add_argument('-p','--multichromosome',type=int,default=0,help= 'Used to determine if the target strain is multichromosome; if the strain is multichromosome, set to 1 and is 0 by default.')
parser.add_argument('-gbc','--multichromosomeGenbank',type=str,default='',help= 'If the strain is multichromosome, enter the Genbank filenames for the chromosomes in order, separated by a "/" symbol.')
args = parser.parse_args()


####基础设置
print('-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())  + '\n' +'Run IFP-PPI' + '\n' )
now_path = os.getcwd()
now_seq = os.sep


if args.multichromosome == 0: #如果不是多染色体
    ####交互频率和线性距离字典信息创建
    # 基因位置信息、基因组长度、基因方向、基因名替换
    print('-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())  + '\n' + 'Convert interaction frequency' + '\n' )
    gene_location, genome_length, gene_direction, gene_locus_tag = genbank_package.genbank_location(args.genbank, now_path, now_seq,args.input)
    # 基因名列表、基因和基因互作的基因名列表、基因对和对应的线性距离
    gene_name, gene_gene_interaction_name, gene_distance = gene_distance_package.gene_gene_distance(genome_length, gene_location)
    # 将交互频率按照从小到大进行排序，基因的长度，基因和对应的bin序号，基因对和对应的交互频率
    gene_length, gene1_bin, gene_interaction = bin_interaction_transiform_gene_package.bin_interaction_transiform_gene(genome_length, args.IFmatrix, args.IFmatrixResolution, gene_location, gene_name, gene_gene_interaction_name, now_path,now_seq,args.binBeginNumber,args.input)



    ####对交互频进行排序
    # 排序交互频率,交互频率从小到大排序
    gene_interaction_sort = sorted(gene_interaction.items(), key=lambda x: x[1])  # 已经是列表了


    ####将排序的基因对按照前10000,排序获取
    # 创建函数，获取基因对字典
    def front_gene_pair_dict(gene_pair_sort, count):
        gene_pair_front = {}
        gene_pair_list = gene_pair_sort[-count:]
        for i in gene_pair_list:
            gene_pair_front[i[0]] = i[1]

        return gene_pair_front
    # 交互频率前10000基因对
    print('-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())  + '\n' + 'Select gene pairs' + '\n' )
    gene_interaction_remove_sort_10000 = front_gene_pair_dict(gene_interaction_sort, args.genepairNumber)


    ####选取前10000的基因对，然后把线性距离小于300bp的进行删除,并且把基因对的数据库信息加进去，也就是在不在string数据库中,最后将基因导出，按照cid和不同的格式
    def filter_final_gene(gene_distance_np, distance_remove, gene_interaction_np, gene_sort_dict):
        """
    将排序前10000的基因对中线性距离小于300的基因对删除
        :param gene_distance: 基因对和线性距离的字典
        :param distance_remove: 删除限制的线性距离
        :param gene_interaction: 交互频率的字典
        :param gene_sort_dict: 前10000基因对字典
        database: 数据库字典
        """
        gene_pair_information = {}
        for i in gene_sort_dict.keys():
            y = []
            if gene_distance_np[i] >= distance_remove:
                y.append(gene_interaction_np[i])
                y.append(gene_distance_np[i])
                gene_pair_information[i] = y
        return gene_pair_information
    #基因对的信息列表是交互频率、线性距离、是否在string数据库中，并且删除了300bp以内的基因对
    print('-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())  + '\n' + 'Remove close gene pairs'  + '\n' )
    gene_pair_remove_10000_information = filter_final_gene(gene_distance,args.removeDistance,gene_interaction,gene_interaction_remove_sort_10000)
    #挑选的基因对中的基因进行提取
    def gene_single(gene_pair_1):
        gene_name = []
        for i in gene_pair_1.keys():
            i1 = i.split(',',1)
            gene_name.append(i1[0])
            gene_name.append(i1[1])
        gene_single_name = list(set(gene_name))
        return  gene_single_name
    gene_pair_remove_10000_information_single = gene_single(gene_pair_remove_10000_information)


    ####如果不识别CIDs
    if args.CIDmatrix == 'None': #如果设置为None，则不进行CID的识别
        now_path_ppi = now_path + now_seq + args.output + now_seq + 'network'
        if not os.path.exists(now_path_ppi):
            os.makedirs(now_path_ppi)
        os.chdir('%s%s%s%snetwork' % (now_path, now_seq, args.output,now_seq))
        with open('gene_association_network.txt', 'w') as file3:
            file3.write('node1\tnode2\tinteraction\tdistance\n')
            for keys3, value3 in gene_pair_remove_10000_information.items():
                keys3_1 = keys3.split(',', -1)
                file3.write(keys3_1[0] + '\t' + keys3_1[1] + '\t' + str(value3[0]) + '\t' + str(int(value3[1])) + '\n')
    else:
        ####创建文件夹
        now_path_cid1 = now_path + now_seq + args.output + now_seq + 'cid' + now_seq + 'cid_range'
        now_path_cid2 = now_path + now_seq + args.output + now_seq + 'cid' + now_seq + 'identify_cid'
        if not os.path.exists(now_path_cid1):
            os.makedirs(now_path_cid1)
        if not os.path.exists(now_path_cid2):
            os.makedirs(now_path_cid2)
        now_path_ppi = now_path + now_seq + args.output + now_seq + 'network'
        if not os.path.exists(now_path_ppi):
            os.makedirs(now_path_ppi)
        if args.CIDmatrix == '': # 如果没有设置CID的矩阵文件名
            # 既没有CIDmatrix，则直接使用IFmatrix; 没有CIDmatrixResolution,则直接使用IFmatrixResolution
            ####识别CID
            # CID的位置信息
            print('-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())  + '\n' + 'Identify CIDs'  + '\n' )
            cid_location, bin_ttest_p_cid, cid_location_number, bin_location = cid_divide_package.bin_interaction_transiform_gene_cid(genome_length, args.IFmatrix, args.IFmatrixResolution, now_path, now_seq, args.binBeginNumber,args.input)
            # 把基因按照CID划分
            cid_gene = cid_total_package.cid_gene(cid_location, gene_location, gene_name, args.IFmatrixResolution)
            # 获取cid的边界和范围
            def cid_scale_amend(cid_location_number_np, bin_location_np):
                cid_scale = {}
                if cid_location_number_np[-1][-1] == list(bin_location_np.keys())[-1]:
                    j = 1
                    for i in range(len(cid_location_number_np)):
                        cid_scale[j] = cid_location_number_np[j]
                        j = j + 1
                else:
                    l = 2
                    cid_scale[1] = [cid_location_number_np[-1][-1] + 1, cid_location_number_np[0][-1]]
                    for i in range(1, len(cid_location_number_np)):
                        cid_scale[l] = cid_location_number_np[i]
                        l = l + 1
                return cid_scale

            cid_scale = cid_scale_amend(cid_location_number, bin_location)
            # 判断基因的cid信息
            gene_pair_remove_10000_information_single_cid = {}  # 选中的基因的cid信息
            for i in cid_gene.keys():
                cid_gene_list = []
                for j in gene_pair_remove_10000_information_single:
                    if j in cid_gene[i]:
                        cid_gene_list.append(j)
                gene_pair_remove_10000_information_single_cid[i] = cid_gene_list
            os.chdir('%s%s%s%scid%scid_range' % (now_path, now_seq, args.output, now_seq, now_seq))
            # 输出cid的范围
            with open('cid_range.txt', 'w') as g1:
                g1.write('CID\tApproximate Genome Position (' + str(args.IFmatrixResolution) + 'bp)\n')
                for i, j in cid_scale.items():
                    g1.write(str(i) + '\t' + str(j[0]) + '-' + str(j[1]) + '\n')
            # 输出cid的边界
            with open('cid_boundary.txt', 'w') as g1:
                g1.write('cid_boundary\n')
                for i in cid_scale.values():
                    g1.write(str(i[1]) + '\n')
            # 把得到的每个bin的t值和p值输出
            os.chdir('%s%s%s%scid%sidentify_cid' % (now_path, now_seq, args.output, now_seq, now_seq))
            with open('bin_t_p.txt', 'w') as p1:
                p1.write('T_value\tP_value\tBin\tDirection\n')
                for p2, p3 in bin_ttest_p_cid.items():
                    if p3[0] > 0:
                        p1.write('%s\t%s\t%s\tLight\n' % (str(p3[0]), str(p3[1]), p2))
                    else:
                        p1.write('%s\t%s\t%s\tRight\n' % (str(p3[0]), str(p3[1]), p2))
            # 输出蛋白互作的基因对，绘制蛋白互作对
            os.chdir('%s%s%s%snetwork' % (now_path, now_seq, args.output, now_seq))
            with open('gene_association_network.txt', 'w') as file3:
                file3.write('node1\tnode2\tinteraction\tdistance\n')
                for keys3, value3 in gene_pair_remove_10000_information.items():
                    keys3_1 = keys3.split(',', -1)
                    file3.write(keys3_1[0] + '\t' + keys3_1[1] + '\t' + str(value3[0]) + '\t' + str(int(value3[1])) + '\n')
            with open('gene_cid.txt', 'w') as file4:
                file4.write('node\tcid\n')
                for i in gene_pair_remove_10000_information_single_cid.keys():
                    for j in gene_pair_remove_10000_information_single_cid[i]:
                        file4.write(j + '\t' + i + '\n')
            with open('gene_all_cid.txt','w')as file5:
                file5.write('gene\tcid\n')
                for i in cid_gene.keys():
                    for j in cid_gene[i]:
                        file5.write(j + '\t' + i + '\n')
        else:
            ####识别CID
            # CID的位置信息
            print('-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())  + '\n' + 'Identify CIDs' +  '\n' )
            cid_location, bin_ttest_p_cid, cid_location_number, bin_location = cid_divide_package.bin_interaction_transiform_gene_cid(genome_length, args.CIDmatrix, args.CIDmatrixResolution, now_path, now_seq, args.binBeginNumber,args.input)
            # 把基因按照CID划分
            cid_gene = cid_total_package.cid_gene(cid_location, gene_location, gene_name, args.CIDmatrixResolution)
            # 获取cid的边界和范围
            def cid_scale_amend(cid_location_number_np, bin_location_np):
                cid_scale = {}
                if cid_location_number_np[-1][-1] == list(bin_location_np.keys())[-1]:
                    j = 1
                    for i in range(len(cid_location_number_np)):
                        cid_scale[j] = cid_location_number_np[j]
                        j = j + 1
                else:
                    l = 2
                    cid_scale[1] = [cid_location_number_np[-1][-1] + 1, cid_location_number_np[0][-1]]
                    for i in range(1, len(cid_location_number_np)):
                        cid_scale[l] = cid_location_number_np[i]
                        l = l + 1
                return cid_scale
            cid_scale = cid_scale_amend(cid_location_number, bin_location)
            # 判断基因的cid信息
            gene_pair_remove_10000_information_single_cid = {}  # 选中的基因的cid信息
            for i in cid_gene.keys():
                cid_gene_list = []
                for j in gene_pair_remove_10000_information_single:
                    if j in cid_gene[i]:
                        cid_gene_list.append(j)
                gene_pair_remove_10000_information_single_cid[i] = cid_gene_list
            os.chdir('%s%s%s%scid%scid_range' % (now_path, now_seq, args.output, now_seq, now_seq))
            # 输出cid的范围
            with open('cid_range.txt', 'w') as g1:
                g1.write('CID\tApproximate Genome Position (' + str(args.CIDmatrixResolution) + 'bp)\n')
                for i, j in cid_scale.items():
                    g1.write(str(i) + '\t' + str(j[0]) + '-' + str(j[1]) + '\n')
            # 输出cid的边界
            with open('cid_boundary.txt', 'w') as g1:
                g1.write('cid_boundary\n')
                for i in cid_scale.values():
                    g1.write(str(i[1]) + '\n')
            # 把得到的每个bin的t值和p值输出
            os.chdir('%s%s%s%scid%sidentify_cid' % (now_path, now_seq, args.output, now_seq, now_seq))
            # 把得到的bin和对应的左右的log2打分对比
            with open('bin_t_p.txt', 'w') as p1:
                p1.write('T_value\tP_value\tBin\tDirection\n')
                for p2, p3 in bin_ttest_p_cid.items():
                    if p3[0] > 0:
                        p1.write('%s\t%s\t%s\tLight\n' % (str(p3[0]), str(p3[1]), p2))
                    else:
                        p1.write('%s\t%s\t%s\tRight\n' % (str(p3[0]), str(p3[1]), p2))
            # 输出蛋白互作的基因对，绘制蛋白互作对
            os.chdir('%s%s%s%snetwork' % (now_path, now_seq, args.output, now_seq))
            with open('gene_association_network.txt', 'w') as file3:
                file3.write('node1\tnode2\tinteraction\tdistance\n')
                for keys3, value3 in gene_pair_remove_10000_information.items():
                    keys3_1 = keys3.split(',', -1)
                    file3.write(keys3_1[0] + '\t' + keys3_1[1] + '\t' + str(value3[0]) + '\t' + str(int(value3[1])) + '\n')
            with open('gene_cid.txt', 'w') as file4:
                file4.write('node\tcid\n')
                for i in gene_pair_remove_10000_information_single_cid.keys():
                    for j in gene_pair_remove_10000_information_single_cid[i]:
                        file4.write(j + '\t' + i + '\n')
            with open('gene_all_cid.txt','w')as file5:
                file5.write('gene\tcid\n')
                for i in cid_gene.keys():
                    for j in cid_gene[i]:
                        file5.write(j + '\t' + i + '\n')
else: #多染色体
    genebank_name = args.multichromosomeGenbank.split('/', -1)  # genbank文件名列表
    gene_location = {}
    genome_length = 0
    gene_direction = {}
    gene_locus_tag = {}
    gene_name = []
    gene_distance = {}
    gene_length = {}
    gene1_bin = {}
    gene_chr = {}  # 不同染色体的基因
    gene_total_bin = []  # 不同染色体的bin序号
    gene_n = {}
    def Merge(dict1, dict2):
        res = {**dict1, **dict2}
        return res
    print('-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S',time.localtime()) + '\n' + 'Convert interaction frequency' + '\n')
    for i in range(len(genebank_name)):
        gene_location_sub, genome_length_sub, gene_direction_sub, gene_locus_tag_sub = genbank_package.genbank_location(genebank_name[i], now_path, now_seq, args.input)
        gene_name_sub, gene_gene_interaction_name_sub, gene_distance_sub = gene_distance_package.gene_gene_distance(genome_length_sub, gene_location_sub)
        if i > 0:  # 不是第一条染色体,其中bin_count是bin的数量,gene1_bin_sub中也是这条染色体上的bin的排序,这里就一直是1
            gene_length_sub, gene1_bin_sub, bin_count, gene_n_sub = bin_interaction_transiform_gene_package.bin_interaction_transiform_gene_poly(
                genome_length_sub, args.IFmatrixResolution, gene_location_sub, gene_name_sub, 1)
        else:  # 第一条染色体，其中是bin_count是第一条染色体的最后一个bin的序号，这是是bin开始的序号
            gene_length_sub, gene1_bin_sub, bin_count, gene_n_sub = bin_interaction_transiform_gene_package.bin_interaction_transiform_gene_poly(
                genome_length_sub, args.IFmatrixResolution, gene_location_sub, gene_name_sub, args.binBeginNumber)
        gene_location = Merge(gene_location, gene_location_sub)
        genome_length = genome_length + genome_length_sub
        gene_direction = Merge(gene_direction, gene_direction_sub)
        gene_locus_tag = Merge(gene_locus_tag, gene_locus_tag_sub)
        gene_name = gene_name + gene_name_sub
        gene_distance = Merge(gene_distance, gene_distance_sub)
        gene_length = Merge(gene_length, gene_length_sub)
        gene1_bin = Merge(gene1_bin, gene1_bin_sub)
        gene_chr[i] = gene_name_sub
        gene_total_bin.append(bin_count)
        gene_n = Merge(gene_n, gene_n_sub)

    # 计算gene1_bin，总的
    gene1_bin_merge = {}
    for j in gene_chr.keys():
        if j != 0:
            for r in gene_chr[j]:  # 遍历所有在第j条染色体的基因
                bin_list = []  # 新的bin的序列列表
                for u in gene1_bin[r]:  # 遍历第j条染色体基因的bin序号
                    for h in range(j):  # 遍历第j条染色体之前的bin序号，gene_total_bin[0]是第1条的最后bin的序号，gene_total_bin[1]是第2条的bin的数量
                        l_bin1 = int(u) + gene_total_bin[h]  # 遍历相加以后就是现在基因对应的bin的序号
                        bin_list.append(str(l_bin1))
                gene1_bin_merge[r] = bin_list
        else:
            for k in gene_chr[0]:  # 如果是在第一条染色体，直接不动，直接转换过来
                gene1_bin_merge[k] = gene1_bin[k]

    gene_gene_interaction_name = []  # 这个需要根据最后总的互作基因名进行获取
    for i_gene_name in range(len(gene_name)):
        i_gene_name_2 = i_gene_name + 1  # 保证比i_gene_name大一个，防止有一样的互作名
        while i_gene_name_2 < len(gene_name):
            i_gene_name_3 = '%s,%s' % (gene_name[i_gene_name], gene_name[i_gene_name_2])
            gene_gene_interaction_name.append(i_gene_name_3)
            i_gene_name_2 = i_gene_name_2 + 1

    # 计算交互频率
    gene_interaction = bin_interaction_transiform_gene_package.bin_interaction_transiform_gene_interaction(gene_n,gene1_bin_merge,args.IFmatrix,gene_gene_interaction_name,now_path,now_seq,args.input)

    ####对交互频进行排序
    # 排序交互频率,交互频率从小到大排序
    gene_interaction_sort = sorted(gene_interaction.items(), key=lambda x: x[1])  # 已经是列表了


    ####将排序的基因对按照前10000,排序获取
    # 创建函数，获取基因对字典
    def front_gene_pair_dict(gene_pair_sort, count):
        gene_pair_front = {}
        gene_pair_list = gene_pair_sort[-count:]
        for i in gene_pair_list:
            gene_pair_front[i[0]] = i[1]

        return gene_pair_front


    # 交互频率前10000基因对
    print('-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()) + '\n' + 'Select gene pairs' + '\n')
    gene_interaction_remove_sort_10000 = front_gene_pair_dict(gene_interaction_sort, args.genepairNumber)


    ####选取前10000的基因对，然后把线性距离小于300bp的进行删除,并且把基因对的数据库信息加进去，也就是在不在string数据库中,最后将基因导出，按照cid和不同的格式
    def filter_final_gene(gene_distance_np, distance_remove, gene_interaction_np, gene_sort_dict):
        """
    将排序前10000的基因对中线性距离小于300的基因对删除
        :param gene_distance: 基因对和线性距离的字典
        :param distance_remove: 删除限制的线性距离
        :param gene_interaction: 交互频率的字典
        :param gene_sort_dict: 前10000基因对字典
        database: 数据库字典
        """
        gene_pair_information = {}
        for i in gene_sort_dict.keys():
            y = []
            if i in gene_distance_np.keys():
                if gene_distance_np[i] >= distance_remove:
                    y.append(gene_interaction_np[i])
                    y.append(gene_distance_np[i])
                    gene_pair_information[i] = y
            else:
                y.append(gene_interaction_np[i])
                y.append('-')
                gene_pair_information[i] = y
        return gene_pair_information


    # 基因对的信息列表是交互频率、线性距离、是否在string数据库中，并且删除了300bp以内的基因对
    print('-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S',
                                          time.localtime()) + '\n' + 'Remove close gene pairs' + '\n')
    gene_pair_remove_10000_information = filter_final_gene(gene_distance, args.removeDistance, gene_interaction,
                                                           gene_interaction_remove_sort_10000)


    # 挑选的基因对中的基因进行提取
    def gene_single(gene_pair_1):
        gene_name = []
        for i in gene_pair_1.keys():
            i1 = i.split(',', 1)
            gene_name.append(i1[0])
            gene_name.append(i1[1])
        gene_single_name = list(set(gene_name))
        return gene_single_name


    gene_pair_remove_10000_information_single = gene_single(gene_pair_remove_10000_information)

    ####如果不识别CIDs
    if args.CIDmatrix == 'None':  # 如果设置为None，则不进行CID的识别
        now_path_ppi = now_path + now_seq + args.output + now_seq + 'network'
        if not os.path.exists(now_path_ppi):
            os.makedirs(now_path_ppi)
        os.chdir('%s%s%s%snetwork' % (now_path, now_seq, args.output, now_seq))
        with open('gene_association_network.txt', 'w') as file3:
            file3.write('node1\tnode2\tinteraction\tdistance\n')
            for keys3, value3 in gene_pair_remove_10000_information.items():
                keys3_1 = keys3.split(',', -1)
                file3.write(keys3_1[0] + '\t' + keys3_1[1] + '\t' + str(value3[0]) + '\t' + str(value3[1]) + '\n')
    else:
        ####创建文件夹
        now_path_cid1 = now_path + now_seq + args.output + now_seq + 'cid' + now_seq + 'cid_range'
        now_path_cid2 = now_path + now_seq + args.output + now_seq + 'cid' + now_seq + 'identify_cid'
        if not os.path.exists(now_path_cid1):
            os.makedirs(now_path_cid1)
        if not os.path.exists(now_path_cid2):
            os.makedirs(now_path_cid2)
        now_path_ppi = now_path + now_seq + args.output + now_seq + 'network'
        if not os.path.exists(now_path_ppi):
            os.makedirs(now_path_ppi)
        if args.CIDmatrix == '':  # 如果没有设置CID的矩阵文件名
            # 既没有CIDmatrix，则直接使用IFmatrix; 没有CIDmatrixResolution,则直接使用IFmatrixResolution
            # CID的位置信息
            print( '-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()) + '\n' + 'Identify CIDs' + '\n')
            genebank_name = args.multichromosomeGenbank.split('/', -1)  # genbank文件名列表
            cid_count = 0  # 这是CID的数量
            # 多条染色体的各自的bin_location
            bin_count_all = 0
            cid_gene_total = {}
            def Merge(dict1, dict2):
                res = {**dict1, **dict2}
                return res
            for chr_number in range(len(genebank_name)):
                gene_location_sub, genome_length_sub, gene_direction_sub, gene_locus_tag_sub = genbank_package.genbank_location(genebank_name[chr_number], now_path, now_seq, args.input)
                gene_name_sub, gene_gene_interaction_name_sub, gene_distance_sub = gene_distance_package.gene_gene_distance(genome_length_sub, gene_location_sub)
                if chr_number == 0:  # 这是是bin的初始序号
                    bin_location, bin_light_right, bin_count, bin_count_begin = cid_divide_package.bin_location_cid(
                        genome_length_sub, args.IFmatrixResolution, args.binBeginNumber)
                    cid_location, bin_ttest_p_cid, cid_location_number = cid_divide_package.gene_cid_divide(
                        bin_location, args.IFmatrix, bin_light_right, now_path, now_seq, args.input, args.binBeginNumber)
                    os.chdir('%s%s%s%scid%sidentify_cid' % (now_path, now_seq, args.output, now_seq, now_seq))
                    # 把得到的bin和对应的左右的log2打分对比
                    with open('bin_t_p_chr' + str(chr_number+1) + '.txt', 'w') as p1:
                        p1.write('T_value\tP_value\tBin\tDirection\n')
                        for p2, p3 in bin_ttest_p_cid.items():
                            if p3[0] > 0:
                                p1.write('%s\t%s\t%s\tLight\n' % (str(p3[0]), str(p3[1]), str(int(p2) - bin_count_all)))
                            else:
                                p1.write('%s\t%s\t%s\tRight\n' % (str(p3[0]), str(p3[1]), str(int(p2) - bin_count_all)))
                else:  # 这全是1.计算的是数量
                    bin_location, bin_light_right, bin_count, bin_count_begin = cid_divide_package.bin_location_cid2(
                        genome_length_sub, args.IFmatrixResolution, bin_count_all)
                    cid_location, bin_ttest_p_cid, cid_location_number = cid_divide_package.gene_cid_divide2(
                        bin_location, args.IFmatrix, bin_light_right, now_path, now_seq, args.input, bin_count_begin,
                        bin_count_all, args.IFmatrixResolution)
                    os.chdir('%s%s%s%scid%sidentify_cid' % (now_path, now_seq, args.output, now_seq, now_seq))
                    # 把得到的bin和对应的左右的log2打分对比
                    with open('bin_t_p_chr' + str(chr_number+1) + '.txt', 'w') as p1:
                        p1.write('T_value\tP_value\tBin\tDirection\n')
                        for p2, p3 in bin_ttest_p_cid.items():
                            if p3[0] > 0:
                                p1.write('%s\t%s\t%s\tLight\n' % (str(p3[0]), str(p3[1]), str(int(p2) - bin_count_all)))
                            else:
                                p1.write('%s\t%s\t%s\tRight\n' % (str(p3[0]), str(p3[1]), str(int(p2) - bin_count_all)))


                # 获取cid的边界和范围
                def cid_scale_amend(cid_location_number_np, bin_location_np):
                    cid_scale = {}
                    if cid_location_number_np[-1][-1] == list(bin_location_np.keys())[-1]:
                        j = 1
                        for i in range(len(cid_location_number_np)):
                            cid_scale[j] = cid_location_number_np[j]
                            j = j + 1
                    else:
                        l = 2
                        cid_scale[1] = [cid_location_number_np[-1][-1] + 1, cid_location_number_np[0][-1]]
                        for i in range(1, len(cid_location_number_np)):
                            cid_scale[l] = cid_location_number_np[i]
                            l = l + 1
                    return cid_scale

                cid_scale = cid_scale_amend(cid_location_number, bin_location)
                cid_gene = cid_total_package.cid_gene2(cid_location, gene_location_sub, gene_name_sub,args.IFmatrixResolution,cid_count,chr_number)
                # 判断基因的cid信息
                os.chdir('%s%s%s%scid%scid_range' % (now_path, now_seq, args.output, now_seq, now_seq))
                # 输出cid的范围
                with open('cid_range_chr' + str(chr_number+1) + '.txt', 'w') as g1:
                    g1.write('CID\tApproximate Genome Position (' + str(args.IFmatrixResolution) + 'bp)\n')
                    for i, j in cid_scale.items():
                        g1.write(str(i + cid_count) + '\t' + str(j[0] - bin_count_all) + '-' + str(j[1] - bin_count_all) + '\n')
                # 输出cid的边界
                with open('cid_boundary_chr' + str(chr_number+1) + '.txt', 'w') as g1:
                    g1.write('cid_boundary\n')
                    for i in cid_scale.values():
                        g1.write(str(i[1] - bin_count_all) + '\n')
                cid_count = cid_count + len(cid_scale)
                bin_count_all = bin_count + bin_count_all
                cid_gene_total = Merge(cid_gene,cid_gene_total)
            # 判断基因的cid信息
            gene_pair_remove_10000_information_single_cid = {}  # 选中的基因的cid信息
            for i in cid_gene_total.keys():
                cid_gene_list = []
                for j in gene_pair_remove_10000_information_single:
                    if j in cid_gene_total[i]:
                        cid_gene_list.append(j)
                gene_pair_remove_10000_information_single_cid[i] = cid_gene_list
            # 输出蛋白互作的基因对，绘制蛋白互作对
            os.chdir('%s%s%s%snetwork' % (now_path, now_seq, args.output, now_seq))
            with open('gene_association_network.txt', 'w') as file3:
                file3.write('node1\tnode2\tinteraction\tdistance\n')
                for keys3, value3 in gene_pair_remove_10000_information.items():
                    keys3_1 = keys3.split(',', -1)
                    file3.write(
                        keys3_1[0] + '\t' + keys3_1[1] + '\t' + str(value3[0]) + '\t' + str(value3[1]) + '\n')
            with open('gene_cid.txt', 'w') as file4:
                file4.write('node\tcid\n')
                for i in gene_pair_remove_10000_information_single_cid.keys():
                    for j in gene_pair_remove_10000_information_single_cid[i]:
                        file4.write(j + '\t' + i + '\n')
            with open('gene_all_cid.txt', 'w') as file5:
                file5.write('gene\tcid\n')
                for i in cid_gene_total.keys():
                    for j in cid_gene_total[i]:
                        file5.write(j + '\t' + i + '\n')
        else:
            ####识别CID
            # CID的位置信息
            print( '-' * 30 + '\n' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()) + '\n' + 'Identify CIDs' + '\n')
            genebank_name = args.multichromosomeGenbank.split('/', -1)  # genbank文件名列表
            cid_count = 0  # 这是CID的数量
            # 多条染色体的各自的bin_location
            bin_count_all = 0
            cid_gene_total = {}
            def Merge(dict1, dict2):
                res = {**dict1, **dict2}
                return res
            for chr_number in range(len(genebank_name)):
                gene_location_sub, genome_length_sub, gene_direction_sub, gene_locus_tag_sub = genbank_package.genbank_location(genebank_name[chr_number], now_path, now_seq, args.input)
                gene_name_sub, gene_gene_interaction_name_sub, gene_distance_sub = gene_distance_package.gene_gene_distance(genome_length_sub, gene_location_sub)
                if chr_number == 0:  # 这是是bin的初始序号
                    bin_location, bin_light_right, bin_count, bin_count_begin = cid_divide_package.bin_location_cid(
                        genome_length_sub, args.CIDmatrixResolution, args.binBeginNumber)
                    cid_location, bin_ttest_p_cid, cid_location_number = cid_divide_package.gene_cid_divide(
                        bin_location, args.CIDmatrix, bin_light_right, now_path, now_seq, args.input, args.binBeginNumber)
                    os.chdir('%s%s%s%scid%sidentify_cid' % (now_path, now_seq, args.output, now_seq, now_seq))
                    # 把得到的bin和对应的左右的log2打分对比
                    with open('bin_t_p_chr' + str(chr_number+1) + '.txt', 'w') as p1:
                        p1.write('T_value\tP_value\tBin\tDirection\n')
                        for p2, p3 in bin_ttest_p_cid.items():
                            if p3[0] > 0:
                                p1.write('%s\t%s\t%s\tLight\n' % (str(p3[0]), str(p3[1]), str(int(p2) - bin_count_all)))
                            else:
                                p1.write('%s\t%s\t%s\tRight\n' % (str(p3[0]), str(p3[1]), str(int(p2) - bin_count_all)))
                else:  # 这全是1.计算的是数量
                    bin_location, bin_light_right, bin_count, bin_count_begin = cid_divide_package.bin_location_cid2(
                        genome_length_sub, args.CIDmatrixResolution, bin_count_all)
                    cid_location, bin_ttest_p_cid, cid_location_number = cid_divide_package.gene_cid_divide2(
                        bin_location, args.CIDmatrix, bin_light_right, now_path, now_seq, args.input, bin_count_begin,
                        bin_count_all, args.CIDmatrixResolution)
                    os.chdir('%s%s%s%scid%sidentify_cid' % (now_path, now_seq, args.output, now_seq, now_seq))
                    # 把得到的bin和对应的左右的log2打分对比
                    with open('bin_t_p_chr' + str(chr_number+1) + '.txt', 'w') as p1:
                        p1.write('T_value\tP_value\tBin\tDirection\n')
                        for p2, p3 in bin_ttest_p_cid.items():
                            if p3[0] > 0:
                                p1.write('%s\t%s\t%s\tLight\n' % (str(p3[0]), str(p3[1]), str(int(p2) - bin_count_all)))
                            else:
                                p1.write('%s\t%s\t%s\tRight\n' % (str(p3[0]), str(p3[1]), str(int(p2) - bin_count_all)))


                # 获取cid的边界和范围
                def cid_scale_amend(cid_location_number_np, bin_location_np):
                    cid_scale = {}
                    if cid_location_number_np[-1][-1] == list(bin_location_np.keys())[-1]:
                        j = 1
                        for i in range(len(cid_location_number_np)):
                            cid_scale[j] = cid_location_number_np[j]
                            j = j + 1
                    else:
                        l = 2
                        cid_scale[1] = [cid_location_number_np[-1][-1] + 1, cid_location_number_np[0][-1]]
                        for i in range(1, len(cid_location_number_np)):
                            cid_scale[l] = cid_location_number_np[i]
                            l = l + 1
                    return cid_scale

                cid_scale = cid_scale_amend(cid_location_number, bin_location)
                cid_gene = cid_total_package.cid_gene2(cid_location, gene_location_sub, gene_name_sub,args.CIDmatrixResolution,cid_count,chr_number)
                # 判断基因的cid信息
                os.chdir('%s%s%s%scid%scid_range' % (now_path, now_seq, args.output, now_seq, now_seq))
                # 输出cid的范围
                with open('cid_range_chr' + str(chr_number + 1) + '.txt', 'w') as g1:
                    g1.write('CID\tApproximate Genome Position (' + str(args.CIDmatrixResolution) + 'bp)\n')
                    for i, j in cid_scale.items():
                        g1.write(str(i + cid_count) + '\t' + str(j[0] - bin_count_all) + '-' + str(
                            j[1] - bin_count_all) + '\n')
                # 输出cid的边界
                with open('cid_boundary_chr' + str(chr_number + 1) + '.txt', 'w') as g1:
                    g1.write('cid_boundary\n')
                    for i in cid_scale.values():
                        g1.write(str(i[1] - bin_count_all) + '\n')
                cid_count = cid_count + len(cid_scale)
                bin_count_all = bin_count + bin_count_all
                cid_gene_total = Merge(cid_gene,cid_gene_total)
            # 判断基因的cid信息
            gene_pair_remove_10000_information_single_cid = {}  # 选中的基因的cid信息
            for i in cid_gene_total.keys():
                cid_gene_list = []
                for j in gene_pair_remove_10000_information_single:
                    if j in cid_gene_total[i]:
                        cid_gene_list.append(j)
                gene_pair_remove_10000_information_single_cid[i] = cid_gene_list
            # 输出蛋白互作的基因对，绘制蛋白互作对
            os.chdir('%s%s%s%snetwork' % (now_path, now_seq, args.output, now_seq))
            with open('gene_association_network.txt', 'w') as file3:
                file3.write('node1\tnode2\tinteraction\tdistance\n')
                for keys3, value3 in gene_pair_remove_10000_information.items():
                    keys3_1 = keys3.split(',', -1)
                    file3.write(
                        keys3_1[0] + '\t' + keys3_1[1] + '\t' + str(value3[0]) + '\t' + str(value3[1]) + '\n')
            with open('gene_cid.txt', 'w') as file4:
                file4.write('node\tcid\n')
                for i in gene_pair_remove_10000_information_single_cid.keys():
                    for j in gene_pair_remove_10000_information_single_cid[i]:
                        file4.write(j + '\t' + i + '\n')
            with open('gene_all_cid.txt', 'w') as file5:
                file5.write('gene\tcid\n')
                for i in cid_gene_total.keys():
                    for j in cid_gene_total[i]:
                        file5.write(j + '\t' + i + '\n')


end = time.time()
print('Time:    ' + str(end - start) + 's') #整个程序运行的时间



