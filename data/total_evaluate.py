####
# 时间：20220220
# 地点：华中农业大学
# 姓名：高远
# 功能：总代码
####


import os
import genbank_evaluate
import gene_distance_evaluate
import bin_interaction_transiform_gene_evaluate
import string_database_evaluate
import prandroc_evaluate
import kegg_evaluate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib



####基础设置
now_path = os.getcwd()
now_seq = os.sep
now_sample_total = ['GSM2870407','GSM2870408','GSM2870409']
matrix_name = ['GSM2870407_1000_iced.matrix','GSM2870408_1000_iced.matrix','GSM2870409_1000_iced.matrix']



####创建文件夹
path1 = ['output']  # 第一层目录
path2 = ['GSM2870407', 'GSM2870408', 'GSM2870409']  # 第二层目录
path3 = ['kegg', 'roc']  # 第三层目录
def mkdir_total(now_path,now_seq,path1,path2,path3): # 把所有的文件夹先创建
    for i in path1:
        for j in path2:
            for l in path3:
                path_new = now_path + now_seq + i + now_seq + j + now_seq + l
                isexist = os.path.exists(path_new)
                if not isexist:
                    os.makedirs(path_new)
                else:
                    n = 'false'
mkdir_total(now_path,now_seq,path1,path2,path3)




for now_sample_list in range(len(now_sample_total)):
    now_sample = now_sample_total[now_sample_list]
    matrix_name_list = matrix_name[now_sample_list]

    ####交互频率和线性距离字典信息创建
    # 基因位置信息、基因组长度、基因方向、基因名替换
    gene_location, genome_length, gene_direction, gene_locus_tag = genbank_evaluate.genbank_location('U00096.3.gb', now_path,now_seq, now_sample)
    # 基因名列表、基因和基因互作的基因名列表、基因对和对应的线性距离
    gene_name, gene_gene_interaction_name, gene_distance = gene_distance_evaluate.gene_gene_distance(genome_length,gene_location)
    # 将交互频率按照从小到大进行排序，基因的长度，基因和对应的bin序号，基因对和对应的交互频率
    gene_length, gene1_bin, gene_interaction = bin_interaction_transiform_gene_evaluate.bin_interaction_transiform_gene(genome_length, matrix_name_list, 1000, gene_location, gene_name, gene_gene_interaction_name,now_path, now_seq, now_sample)


    ####按条件处理交互频率和线性距离字典
    # 删除线性距离中交互频率没有的基因对
    def remove_distance(gene_interaction_np, gene_distance_np):
        gene_distance_new = {}
        for i in gene_interaction_np.keys():
            gene_distance_new[i] = gene_distance_np[i]
        return gene_distance_new


    # 删除线性距离中总的交互频率没有的基因对,删除总的交互频率没有的基因对
    gene_distance_remove = remove_distance(gene_interaction, gene_distance)

    ####对交互频率和线性距离进行排序，然后按照百分比进行划分
    # 排序交互频率和线性距离，交互频率从小到大排序；线性距离从大到小排序；全部基因对
    gene_interaction_sort = sorted(gene_interaction.items(), key=lambda x: x[1])  # 已经是列表了
    gene_distance_remove_sort = sorted(gene_distance_remove.items(), key=lambda x: x[1], reverse=True)


    ####将排序的基因对按照前20000，前10000，前5000排序获取
    # 创建函数，获取基因对字典
    def front_gene_pair_dict(gene_pair_sort, count):
        gene_pair_front = {}
        gene_pair_list = gene_pair_sort[-count:]
        for i in gene_pair_list:
            gene_pair_front[i[0]] = i[1]

        return gene_pair_front


    # 线性距离和交互频率前10000基因对
    gene_distance_remove_sort_10000 = front_gene_pair_dict(gene_distance_remove_sort, 10000)
    gene_interaction_remove_sort_10000 = front_gene_pair_dict(gene_interaction_sort, 10000)

    ####获取数据库数据
    # 基因对在string数据库中得分，string数据库得分字典
    string_dict_total = string_database_evaluate.string_compare(now_path, now_seq, 'K12_MG1655_detailed_v11.5_new.txt',
                                                       now_sample)

    ####计算roc和pr的auc值；string数据库总字典
    # 基因对的标签值，基因对的value值；前10000、线性距离
    gene_true_distance_remove_10000, gene_list_distance_remove_10000 = string_database_evaluate.string_au_distance(now_path,
                                                                                                          now_seq,
                                                                                                          'K12_MG1655_detailed_v11.5_new.txt',
                                                                                                          gene_distance_remove_sort_10000,
                                                                                                          now_sample,
                                                                                                          400)
    # 基因对的标签值，基因对的value值；前10000、交互频率
    gene_true_interaction_remove_10000, gene_list_interaction_remove_10000 = string_database_evaluate.string_au(now_path,
                                                                                                       now_seq,
                                                                                                       'K12_MG1655_detailed_v11.5_new.txt',
                                                                                                       gene_interaction_remove_sort_10000,
                                                                                                       now_sample, 400)
    # 计算aupr和auroc的值;前10000、线性距离
    auroc_distance_remove_10000, fpr_distance_remove_10000, tpr_distance_remove_10000 = prandroc_evaluate.auroc(
        gene_list_distance_remove_10000, gene_true_distance_remove_10000)
    aupr_distance_remove_10000, recall_distance_remove_10000, precision_distance_remove_10000 = prandroc_evaluate.aupr(
        gene_list_distance_remove_10000, gene_true_distance_remove_10000)
    # 计算aupr和auroc的值;前10000、交互频率
    auroc_interaction_remove_10000, fpr_interaction_remove_10000, tpr_interaction_remove_10000 = prandroc_evaluate.auroc(
        gene_list_interaction_remove_10000, gene_true_interaction_remove_10000)
    aupr_interaction_remove_10000, recall_interaction_remove_10000, precision_interaction_remove_10000 = prandroc_evaluate.aupr(
        gene_list_interaction_remove_10000, gene_true_interaction_remove_10000)

    ####计算roc和pr的auc值；KEGG通路
    # kegg数据，通路名和基因名字典
    kegg_dict = kegg_evaluate.kegg_dict(now_path, now_seq, now_sample, 'KEGGREST_WithGene.tsv')


    # 获取所有kegg基因的字典
    def kegg_gene_single(kegg_dict):
        kegg_gene = []
        for i in kegg_dict.values():
            for j in i:
                kegg_gene.append(j)

        kegg_gene_single = list(set(kegg_gene))

        return kegg_gene_single


    kegg_gene = kegg_gene_single(kegg_dict)
    # 制作一个基因对和共同出现通路数的字典,在一个通路里的就算一个，出现在不同通路里的算0个
    gene_kegg_number = {}
    for i in gene_distance.keys():
        n = 0  # 通路数
        i1 = i.split(',', 1)
        for o in kegg_dict.values():
            if i1[0] in o and i1[1] in o:
                n = n + 1
        gene_kegg_number[i] = n


    # 创建函数，用来创建标签和value值，输入基因对字典
    def kegg_number_distance(gene_distance, gene_kegg_number, cutoff):
        gene_true = []
        value_list = []
        for i, i1 in gene_distance.items():
            if i in gene_kegg_number.keys():
                if gene_kegg_number[i] >= cutoff:
                    gene_true.append(1)
                    value_list.append(2320813 - i1)
                else:
                    gene_true.append(0)
                    value_list.append(2320813 - i1)

        return gene_true, value_list


    def kegg_number_interaction(gene_interaction, gene_kegg_number, cutoff):
        gene_true = []
        value_list = []
        for i, i1 in gene_interaction.items():
            if i in gene_kegg_number.keys():
                if gene_kegg_number[i] >= cutoff:
                    gene_true.append(1)
                    value_list.append(i1)
                else:
                    gene_true.append(0)
                    value_list.append(i1)
        return gene_true, value_list


    # 标签值、value值；全部基因对前10000,线性距离
    gene_true_distance_remove_kegg_10000, gene_list_distance_remove_kegg_10000 = kegg_number_distance(
        gene_distance_remove_sort_10000, gene_kegg_number, 1)
    # 标签值、value值；全部基因对前10000,交互频率
    gene_true_interaction_kegg_10000, gene_list_interaction_kegg_10000 = kegg_number_interaction(
        gene_interaction_remove_sort_10000, gene_kegg_number, 1)
    # 计算aupr和auroc的值;全部基因对前10000、线性距离
    auroc_distance_remove_kegg_10000, fpr_distance_remove_kegg_10000, tpr_distance_remove_kegg_10000 = prandroc_evaluate.auroc(
        gene_list_distance_remove_kegg_10000, gene_true_distance_remove_kegg_10000)
    aupr_distance_remove_kegg_10000, recall_distance_remove_kegg_10000, precision_distance_remove_kegg_10000 = prandroc_evaluate.aupr(
        gene_list_distance_remove_kegg_10000, gene_true_distance_remove_kegg_10000)
    # 计算aupr和auroc的值;全部基因对前10000、交互频率
    auroc_interaction_kegg_10000, fpr_interaction_kegg_10000, tpr_interaction_kegg_10000 = prandroc_evaluate.auroc(
        gene_list_interaction_kegg_10000, gene_true_interaction_kegg_10000)
    aupr_interaction_kegg_10000, recall_interaction_kegg_10000, precision_interaction_kegg_10000 = prandroc_evaluate.aupr(
        gene_list_interaction_kegg_10000, gene_true_interaction_kegg_10000)

    '''
    ####画ROC曲线图
    #kegg，前10000对基因对的图
    plt.plot(fpr_distance_remove_kegg_10000,tpr_distance_remove_kegg_10000,color = '#B71C1C',label = 'Distance ROC(area = {0:.3f})'.format(auroc_distance_remove_kegg_10000))
    plt.plot(fpr_interaction_kegg_10000,tpr_interaction_kegg_10000,color = '#00438E',label = 'Interaction ROC(area = {0:.3f})'.format(auroc_interaction_kegg_10000))
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([-0.05, 1.05])  # 设置x、y轴的上下限，以免和边缘重合，更好的观察图像的整体
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')  # 可以使用中文，但需要导入一些库即字体
    plt.title('ROC Curve ')
    plt.legend(loc="lower right")
    os.chdir('%s%soutput%s%s%sroc' % (now_path,now_seq,now_seq,now_sample,now_seq))
    plt.savefig('10000_kegg_roc.pdf',format = 'pdf')
    plt.cla() # 清除axes，即当前 figure 中的活动的axes，但其他axes保持不变。
    plt.clf() # 清除当前 figure 的所有axes，但是不关闭这个 window，所以能继续复用于其他的 plot。

    #string数据库，大于100，前10000对基因对的图
    plt.plot(fpr_distance_remove_10000,tpr_distance_remove_10000,color = '#B71C1C',label = 'Distance ROC(area = {0:.3f})'.format(auroc_distance_remove_10000))
    plt.plot(fpr_interaction_remove_10000,tpr_interaction_remove_10000,color = '#00438E',label = 'Interaction ROC(area = {0:.3f})'.format(auroc_interaction_remove_10000))
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([-0.05, 1.05])  # 设置x、y轴的上下限，以免和边缘重合，更好的观察图像的整体
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')  # 可以使用中文，但需要导入一些库即字体
    plt.title('ROC Curve ')
    plt.legend(loc="lower right")
    os.chdir('%s%soutput%s%s%sroc' % (now_path,now_seq,now_seq,now_sample,now_seq))
    plt.savefig('10000_roc_string.pdf',format = 'pdf')
    plt.cla() # 清除axes，即当前 figure 中的活动的axes，但其他axes保持不变。
    plt.clf() # 清除当前 figure 的所有axes，但是不关闭这个 window，所以能继续复用于其他的 plot。
    ####
    '''

    ####计算10000基因对的线性距离和交互频率的阳性基因对数量
    number_true_genepairs = {}


    def gene_pairs_number(gene_remove_sort_np, string_dict_total_np, gene_distance_np):
        gene_list = []
        for i, y in gene_remove_sort_np.items():
            if i in string_dict_total_np.keys():
                if int(string_dict_total_np[i]) >= 400:
                    gene_list.append(gene_distance_np[i])
        return gene_list


    for number_true_list in [10000]:
        list2 = []
        gene_distance_remove_sort_np = front_gene_pair_dict(gene_distance_remove_sort, number_true_list)
        gene_interaction_remove_sort_np = front_gene_pair_dict(gene_interaction_sort, number_true_list)
        gene_interaction_list = gene_pairs_number(gene_interaction_remove_sort_np, string_dict_total, gene_distance)
        gene_distance_list = gene_pairs_number(gene_distance_remove_sort_np, string_dict_total, gene_distance)
        list2.append(gene_distance_list)
        list2.append(gene_interaction_list)
        number_true_genepairs[number_true_list] = list2
    number_true = {}
    for number_true_list in [10000]:
        list2 = []
        gene_distance_remove_sort_np = front_gene_pair_dict(gene_distance_remove_sort, number_true_list)
        gene_interaction_remove_sort_np = front_gene_pair_dict(gene_interaction_sort, number_true_list)
        gene_true_distance_remove_np, gene_list_distance_remove_np = string_database_package.string_au_distance(now_path,now_seq,'K12_MG1655_detailed_v11.5_new.txt',gene_distance_remove_sort_np,now_sample, 400)
        gene_true_interaction_remove_np, gene_list_interaction_remove_np = string_database_package.string_au(now_path, now_seq,'K12_MG1655_detailed_v11.5_new.txt',gene_interaction_remove_sort_np,now_sample, 400)
        list2.append(gene_true_distance_remove_np.count(1))
        list2.append(gene_true_interaction_remove_np.count(1))
        number_true[number_true_list] = list2
    number_distance = {}
    number_true1 = {}
    # 前10000对
    distance_range2 = [[0], [0, 500], [500, 1000], [1000, 1500], [1500, 2000], [2000, 3000000]]
    for number_distance_list in [10000]:
        list_total = []
        list_distance = []
        list_interaction = []
        list_true = []
        list_range = []
        for distance_range1_list in range(len(distance_range2)):
            number_0_d = 0
            number_0_i = 0
            if distance_range1_list == 0:
                for number_true_genepairs_values1 in number_true_genepairs[number_distance_list][0]:
                    if number_true_genepairs_values1 == 0:
                        number_0_d = number_0_d + 1
                for number_true_genepairs_values2 in number_true_genepairs[number_distance_list][1]:
                    if number_true_genepairs_values2 == 0:
                        number_0_i = number_0_i + 1
            else:
                for number_true_genepairs_values3 in number_true_genepairs[number_distance_list][0]:
                    if distance_range2[distance_range1_list][0] < number_true_genepairs_values3 <= \
                            distance_range2[distance_range1_list][1]:
                        number_0_d = number_0_d + 1
                for number_true_genepairs_values4 in number_true_genepairs[number_distance_list][1]:
                    if distance_range2[distance_range1_list][0] < number_true_genepairs_values4 <= \
                            distance_range2[distance_range1_list][1]:
                        number_0_i = number_0_i + 1
            list_range.append(distance_range2[distance_range1_list])
            list_distance.append(number_0_d)
            list_interaction.append(number_0_i)
        list_total.append(list_distance)
        list_total.append(list_interaction)
        list_total.append(list_range)
        list_true.append(sum(list_distance))
        list_true.append(sum(list_interaction))
        number_distance[number_distance_list] = list_total
        number_true1[number_distance_list] = list_true
    ####
    # kegg的阳性基因数量
    number_true_genepairs_kegg = {}


    def gene_pairs_number_kegg(gene_remove_sort_np, string_dict_total_np, gene_distance_np):
        gene_list = []
        for i, y in gene_remove_sort_np.items():
            if i in string_dict_total_np.keys():
                if int(string_dict_total_np[i]) >= 1:
                    gene_list.append(gene_distance_np[i])
        return gene_list


    for number_true_list in [10000]:
        list2 = []
        gene_distance_remove_sort_np = front_gene_pair_dict(gene_distance_remove_sort, number_true_list)
        gene_interaction_remove_sort_np = front_gene_pair_dict(gene_interaction_sort, number_true_list)
        gene_interaction_list = gene_pairs_number_kegg(gene_interaction_remove_sort_np, gene_kegg_number, gene_distance)
        gene_distance_list = gene_pairs_number_kegg(gene_distance_remove_sort_np, gene_kegg_number, gene_distance)
        list2.append(gene_distance_list)
        list2.append(gene_interaction_list)
        number_true_genepairs_kegg[number_true_list] = list2
    number_distance_kegg = {}
    number_true1_kegg = {}
    # 前10000对
    distance_range2 = [[0], [0, 500], [500, 1000], [1000, 1500], [1500, 2000], [2000, 3000000]]
    for number_distance_list in [10000]:
        list_total = []
        list_distance = []
        list_interaction = []
        list_true = []
        list_range = []
        for distance_range1_list in range(len(distance_range2)):
            number_0_d = 0
            number_0_i = 0
            if distance_range1_list == 0:
                for number_true_genepairs_values1 in number_true_genepairs_kegg[number_distance_list][0]:
                    if number_true_genepairs_values1 == 0:
                        number_0_d = number_0_d + 1
                for number_true_genepairs_values2 in number_true_genepairs_kegg[number_distance_list][1]:
                    if number_true_genepairs_values2 == 0:
                        number_0_i = number_0_i + 1
            else:
                for number_true_genepairs_values3 in number_true_genepairs_kegg[number_distance_list][0]:
                    if distance_range2[distance_range1_list][0] < number_true_genepairs_values3 <= \
                            distance_range2[distance_range1_list][1]:
                        number_0_d = number_0_d + 1
                for number_true_genepairs_values4 in number_true_genepairs_kegg[number_distance_list][1]:
                    if distance_range2[distance_range1_list][0] < number_true_genepairs_values4 <= \
                            distance_range2[distance_range1_list][1]:
                        number_0_i = number_0_i + 1
            list_range.append(distance_range2[distance_range1_list])
            list_distance.append(number_0_d)
            list_interaction.append(number_0_i)
        list_total.append(list_distance)
        list_total.append(list_interaction)
        list_total.append(list_range)
        list_true.append(sum(list_distance))
        list_true.append(sum(list_interaction))
        number_distance_kegg[number_distance_list] = list_total
        number_true1_kegg[number_distance_list] = list_true

    ####输出文件
    # 输出前1000-10000的基因对的阳性基因对数量的线性距离值
    os.chdir('%s%soutput%s%s%sroc' % (now_path, now_seq, now_seq, now_sample, now_seq))
    with open('number_distance_range.txt', 'w') as number_distance_range_file:
        number_distance_range_file.write('Number\tRange\tGroup\tGene_pairs_number\n')
        for number_distance_keys, number_distance_values in number_distance.items():
            for range_number_values in range(len(number_distance_values[0])):
                number_distance_range_file.write(str(number_distance_values[0][range_number_values]) + '\t' + str(
                    number_distance_values[2][range_number_values]) +
                                                 '\t' + 'Distance' + '\t' + str(number_distance_keys) + '\n' + str(
                    number_distance_values[1][range_number_values]) + '\t' + str(
                    number_distance_values[2][range_number_values]) +
                                                 '\t' + 'Interaction' + '\t' + str(number_distance_keys) + '\n')
    # 输出前1000-10000的基因对的阳性基因对数量的线性距离值,kegg
    os.chdir('%s%soutput%s%s%skegg' % (now_path, now_seq, now_seq, now_sample, now_seq))
    with open('number_distance_range.txt', 'w') as number_distance_range_file:
        number_distance_range_file.write('Number\tRange\tGroup\tGene_pairs_number\n')
        for number_distance_keys, number_distance_values in number_distance_kegg.items():
            for range_number_values in range(len(number_distance_values[0])):
                number_distance_range_file.write(str(number_distance_values[0][range_number_values]) + '\t' + str(
                    number_distance_values[2][range_number_values]) +
                                                 '\t' + 'Distance' + '\t' + str(number_distance_keys) + '\n' + str(
                    number_distance_values[1][range_number_values]) + '\t' + str(
                    number_distance_values[2][range_number_values]) +
                                                 '\t' + 'Interaction' + '\t' + str(number_distance_keys) + '\n')

