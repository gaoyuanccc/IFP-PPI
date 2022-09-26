def bin_interaction_transiform_gene( genome_length,interaction_matrix_filename, resolution_fp,gene_location_np1,gene_name_np,gene_gene_interaction_name_np,now_path_np,now_seq_np,now_sample_np):
    """
将bin之间的互作转换为基因间的互作
    :param gene_location_np1: 基因的位置信息的字典
    :param gene_name_np: 基因的名字的形参
    :param gene_gene_interaction_name_np: 互作基因名的列表形参
    :param interaction_bed_filename: Hic-Pro处理得到的数据的bed文件
    :param interaction_matrix_filename: HiC-Pro处理得到的标准化的bin之间的互作频率文件
    :param resolution_fp: HiC-Pro处理是的分辨率
    :return: gene_interaction_sort: 基因间的互作频率，从小到大排列；     gene_length: 基因和其对应的长度；    gene1_bin: 基因和对应的bin序号；     gene_interaction: 基因对和对应的互作频率
    """
    import numpy
    import os
    import pandas as pd
    import numpy as np
    ####导入备注文件，得到bin和它的位置信息的字典
    os.chdir('%s%sinput%s%s%sinteraction_distance' % (now_path_np, now_seq_np, now_seq_np,now_sample_np,now_seq_np))
    bin_count_bed = (genome_length // resolution_fp) + 1
    bin_name_bed = []
    for bin_count_bed_list in range(bin_count_bed):
        bin_name_bed.append(bin_count_bed_list)
    bin_location = {}
    bin_bed_begin = 0
    bin_bed_finna = resolution_fp
    for bin_name_bed_list in bin_name_bed:
        bin_location[str(bin_name_bed_list)] = [bin_bed_begin, bin_bed_finna]
        bin_bed_begin = bin_bed_begin + resolution_fp
        bin_bed_finna = bin_bed_finna + resolution_fp


    ####导入互作数据，得到互作bin序号和互作频率的字典
    data = pd.read_csv(interaction_matrix_filename, sep='\t', header=None)
    bin_interaction = {}
    if data.shape[1] > 3:
        f1 = open(interaction_matrix_filename)
        f2 = f1.readlines()
        data_list = []
        for i in f2:
            i1 = i.strip('\n').strip('\t').split('\t', -1)
            data_list.append(i1)
        data_matrix = pd.DataFrame(np.array(data_list))
        m = data_matrix.shape[0]
        bin_number = []
        for d in range(m):
            bin_number.append(d)
        bin_bin_name = []
        for m_list in range(len(bin_number)):
            m_list2 = m_list
            while m_list2 < len(bin_number):
                m_list3 = '%s,%s' % (bin_number[m_list], bin_number[m_list2])
                bin_bin_name.append(m_list3)
                m_list2 = m_list2 + 1
        for bin_bin_name_list in bin_bin_name:
            bin_bin_name_list1 = bin_bin_name_list.split(',', -1)
            bin_interaction[bin_bin_name_list] = data_matrix.iloc[int(bin_bin_name_list1[0]), int(bin_bin_name_list1[1])]
    else:
        data = pd.read_csv(interaction_matrix_filename, sep='\t', header=None)
        sub = ['bin1', 'bin2', 'interaction']  # 赋予列的名字
        data.columns = sub
        # df = data.set_index(['bin1','bin2'])
        # 转为字典
        bin_interaction = dict(
            [(str(i) + ',' + str(a), str(b)) for i, a, b in zip(data['bin1'], data['bin2'], data['interaction'])])

    ####制作基因的长度字典 字典格式：基因名：基因长度
    gene_length = {}
    for gene_location_key in gene_location_np1.keys():  # 遍历基因名
        y2 = float(gene_location_np1[gene_location_key][1]) - float(
            gene_location_np1[gene_location_key][0])  # 基因的长度对于相应的基因的后面的位置信息减去前面的位置信息
        gene_length[gene_location_key] = y2

    ####制作基因和其对应的bin的信息的字典 字典格式：基因名：bin序号列表
    gene_bin = {}
    for b1 in gene_location_np1.keys():  # 首先遍历基因的名字
        b2 = int(gene_location_np1[b1][0]) // resolution_fp  # 基因的初始位置除于1000，整数部分就是基因的初始位置所在的bin序号
        b3 = int(gene_location_np1[b1][1]) // resolution_fp  # 基因的末尾位置除于1000，整数部分就是基因的末尾位置所在的bin序号
        b4 = []  # 通过添加每个基因相对应的序号列表
        if b2 == b3:  # 当基因的初始位置和末尾位置都位于同一个bin中时
            b4.append(str(b2))  # 那么基因对应的bin就是一个
        else:  # 当基因的初始位置和末尾位置不在同一个bin时
            b5 = int(gene_location_np1[b1][
                         1]) % resolution_fp  # 这是判断是否基因的末尾是否和bin的末尾重合了，如果重合的话，我们还是只算它在一个bin中，所以我们求的是基因的末尾位置比上1000，余数等于零就是重合了
            if b5 == 0:  # 余数等于零的话，说明正好重合了
                while b2 <= b3 - 1:  # 重合的话，我们只算到b3的前一个，那就是b3-1
                    b4.append(str(b2))  # 这时我们只算基因的初始位置所在的bin
                    b2 = b2 + 1
            else:
                while b2 <= b3:  # 我们依次把初始位置位于的bin序号加到末尾位置位于的bin序号就可以了
                    b4.append(str(b2))
                    b2 = b2 + 1
        gene_bin[b1] = b4  # 基因名和相对应的bin序号

    ####得到基因的性质(是一个bin还是多个bin，字典格式 gene: 1或2)和得到基因和其相对应的bin序号列表（字典格式 gene : bin序号列表)
    gene_n = {}
    gene1_bin = {}
    for gene1_name in gene_name_np:  # 获得每个基因的性质和bin序号列表的循环，得到基因的性质(是一个bin还是多个bin，字典格式 gene: 1或2)和得到基因和其相对应的bin序号列表（字典格式 gene : bin序号列表)
        # gene_name：是基因名的列表
        # bin_intraction: 是互作的bin序号和对应的互作频率字典
        # bin_location: 是bin序号对应的位置信息字典
        # gene_bin : 是基因和其对应的所位于的bin序号列表字典
        # gene_length: 是基因和对应的长度字典
        # gene_location: 是基因和其对应的位置信息列表字典
        # gene_n: 是输出的字典，是一个空字典，是基因和其对应的n值
        # gene1_bin: 是输出的字典，是一个空字典，是基因和其对应的有效bin序号列表
        n = 0
        bin_number1 = []
        if len(gene_bin[gene1_name]) == 1:  # 判断基因是否只位于一个bin内
            for a5 in gene_bin[gene1_name]:
                bin_number1.append(a5)  # 基因所位于的bin的列表
            n = 1  # 给基因进行定性，只位于一个bin中
        else:
            if len(gene_bin[gene1_name]) == 2:  # 判断基因是否只位于两个bin内
                o = float(gene_location_np1[gene1_name][1]) - float(
                    bin_location[gene_bin[gene1_name][1]][0])  # o是基因所位于的第二个bin的部分的长度，是由基因的末尾位置减去基因所位于的第二个bin的初始位置
                o1 = gene_length[gene1_name] / o  # o1是基因的长度对o的倍数，是由基因的总长度比上o
                if o1 == 2:  # 判断基因是不是被两个bin平分了
                    for a6 in gene_bin[gene1_name]:
                        bin_number1.append(a6)
                    n = 2  # 说明基因不止位于一个bin内
                else:
                    if o1 > 2:  # 判断基因大于50%的部分
                        bin_number1.append(gene_bin[gene1_name][1])  # 如果大于2说明大于50%的部分在第二个bin
                        n = 1
                    else:
                        bin_number1.append(gene_bin[gene1_name][0])  # 如果小于2说明大于50%的部分在第一个bin
                        n = 1
            else:  # 说明基因所位于的大于等于三个bin
                if gene_location_np1[gene1_name][0] != bin_location[gene_bin[gene1_name][0]][0]:  # 基因的初始位置和第一个bin的初始位置不相同
                    if gene_location_np1[gene1_name][1] != bin_location[gene_bin[gene1_name][-1]][
                        1]:  # 判断基因的末尾位置和最后一个bin的末尾位置不相同
                        for a1 in gene_bin[gene1_name]:  # 如果基因的初始位置和末尾位置都不相同，宾序号就是bin列表去除首尾
                            bin_number1.append(a1)
                        bin_number1.pop(0)
                        bin_number1.pop()
                        if len(bin_number1) == 1:
                            n = 1
                        else:
                            n = 2
                    else:
                        for a2 in gene_bin[gene1_name]:  # bin的序号就是bin序号列表去除第一个
                            bin_number1.append(a2)
                        bin_number1.pop(0)
                        n = 2
                else:
                    if gene_location_np1[gene1_name][1] != bin_location[gene_bin[gene1_name][-1]][
                        1]:  # 判断基因的末尾位置和最后一个bin的末尾位置不同
                        for a3 in gene_bin[gene1_name]:  # bin的序号等于bin序号列表去除最后一个
                            bin_number1.append(a3)
                        bin_number1.pop()
                        n = 2
                    else:
                        for a4 in gene_bin[gene1_name]:  # 那么bin的序号列表就是基因所占得所有bin
                            bin_number1.append(a4)
                        n = 2
        gene_n[gene1_name] = n
        gene1_bin[gene1_name] = bin_number1
    from decimal import Decimal
    ####得到有互作信息基因对和其对应的互作频率值
    gene_remove = 0  # 计算删除的基因对(没有互作信息的)的个数，有校准的作用
    gene_interaction = {}  # 是基因对和其对应的互作频率值
    for d1 in gene_gene_interaction_name_np:
        # gene_gene_interaction_name: 一个互作基因名的列表(不包括自己和自己互作的基因对)列表
        # gene_n : 是基因和其对应的一个有效bin和多个有效bin的字典(1是一个bin，2是多个bin)
        # gene1_bin :是基因和其对应的有效bin的序号列表的字典
        # bin_interaction: 是互作的bin序号和对应的互作频率字典，用来判断有没有互作信息
        # gene_remove : 输出文件，是删除的没有互作信息的基因对的个数
        # gene_interaction : 输出文件，是有互作信息基因对和其对应的互作频率值
        (d2, d3) = d1.split(',', 1)  # 得到互作的基因对名称，gene1是d2,gene2是d3
        gene1_n = gene_n[d2]  # 这是得到gene1是一个bin还是多个bin
        gene2_m = gene_n[d3]  # 这是得到gene2是一个bin还是多个bin
        if gene1_n == 1 and gene2_m == 1:  # 当gene1和gene2都是一个bin的时候
            bin_bin1 = gene1_bin[d2][0] + ',' + gene1_bin[d3][0]  # 这是gene1和gene2所对应的bin的组合
            bin_bin2 = gene1_bin[d3][0] + ',' + gene1_bin[d2][0]  # 这是防止有反着的bin对没被找到
            if bin_bin1 in bin_interaction.keys() or bin_bin2 in bin_interaction.keys():  # 当互作的bin有互作信息时
                if bin_bin1 in bin_interaction.keys():
                    gene_interaction[d1] = float(bin_interaction[bin_bin1])  # 存入字典中，基因1和基因2的互作频率,主要是不知道是正着的有互作信息还是反着的有互作信息
                else:
                    gene_interaction[d1] = float(bin_interaction[bin_bin2])
            else:
                gene_remove = gene_remove + 1  # 不存在的话，去除的基因对加一
        else:  # 当gene1和gene2至少有一个不是一个bin的时候，就是一对多
            bin1_bin2_interaction = []  # gene1和gene2所有有互作信息的bin组合对应的互作频率
            for d4 in gene1_bin[d2]:  # 遍历gene1的bin序号
                for d5 in gene1_bin[d3]:  # 遍历gene2的bin序号
                    d6 = str(d4) + ',' + str(d5)  # 然后把bin组合起来
                    d7 = str(d5) + ',' + str(d4)  # 同样的，防止反着的bin互作没找到
                    if d6 in bin_interaction.keys() or d7 in bin_interaction.keys():  # 当互作的bin有互作信息时
                        if d6 in bin_interaction.keys():  # 存入bin组合对应的互作频率列表中,主要是不知道是正着的有互作信息还是反着的有互作信息
                            bin1_bin2_interaction.append(float(bin_interaction[d6]))  # 存入的是互作频率
                        else:
                            bin1_bin2_interaction.append(float(bin_interaction[d7]))
            if len(bin1_bin2_interaction) == 0:  # 当有互作信息的互作bin个数是0时
                gene_remove = gene_remove + 1  # 删除基因对个数加一
            else:
                gene_interaction[d1] = numpy.mean(bin1_bin2_interaction)  # 存入字典中，基因对对应的是有互作信息的互作频率的平均值
    gene_interaction_decimal = {}
    for g_i_key,g_i_value in gene_interaction.items():
        interaction_d = '%.6f' % g_i_value
        gene_interaction_decimal[g_i_key] = float(interaction_d)
        ####把基因对按照互作频率的大小从小到大进行排序,保存为列表

    return gene_length, gene1_bin, gene_interaction_decimal


