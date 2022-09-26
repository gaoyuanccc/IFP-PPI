####
# 时间：20220305
# 地点：华中农业大学
# 姓名：高远
# 功能：使用GSM2870409的bin10000数据来划分CID
####
def bin_interaction_transiform_gene_cid(genome_length,interaction_matrix_filename, resolution_fp,now_path_np,now_seq_np,bin_begin_number,input_folder):
    """
    将bin之间的互作转换为基因间的互作
    :param gene_location_np1: 基因的位置信息的字典
    :param gene_name_np: 基因的名字的形参
    :param gene_gene_interaction_name_np: 互作基因名的列表形参
    :param interaction_matrix_filename: 标准化的bin之间的互作频率文件
    :param resolution_fp: 识别cid的分辨率
    bin_begin_number: 初始bin的序号
    :return: gene_interaction_sort: 基因间的互作频率，从小到大排列；     gene_length: 基因和其对应的长度；    gene1_bin: 基因和对应的bin序号；     gene_interaction: 基因对和对应的互作频率
    """
    import os
    import pandas as pd
    import numpy as np
    from scipy import stats
    ####导入备注文件，得到bin和它的位置信息的字典
    os.chdir(now_path_np + now_seq_np + input_folder)
    if (genome_length % resolution_fp) == 0:  # 如果没有余值，则不需要除以的商加1
        bin_count_bed = (genome_length // resolution_fp)
    else:
        bin_count_bed = (genome_length // resolution_fp) + 1
    bin_name_bed = []
    for bin_count_bed_list in range(bin_begin_number, bin_count_bed + bin_begin_number):  # bin_number_begin的初始值是0
        bin_name_bed.append(bin_count_bed_list)

    bin_cid_location = {}
    bin_bed_begin = 0
    bin_bed_finna = resolution_fp
    for bin_name_bed_list in bin_name_bed:
        bin_cid_location[str(bin_name_bed_list)] = [bin_bed_begin, bin_bed_finna]
        bin_bed_begin = bin_bed_begin + resolution_fp
        bin_bed_finna = bin_bed_finna + resolution_fp

    ####导入互作数据，得到互作bin序号和互作频率的字典
    data = pd.read_csv(interaction_matrix_filename, sep='\t', header=None)
    bin_cid_interaction = {}
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
            bin_cid_interaction[bin_bin_name_list] = data_matrix.iloc[
                int(bin_bin_name_list1[0]), int(bin_bin_name_list1[1])]
    else:
        data = pd.read_csv(interaction_matrix_filename, sep='\t', header=None)
        sub = ['bin1', 'bin2', 'interaction']  # 赋予列的名字
        data.columns = sub
        # df = data.set_index(['bin1','bin2'])
        # 转为字典
        bin_cid_interaction = dict(
            [(str(i) + ',' + str(a), str(b)) for i, a, b in zip(data['bin1'], data['bin2'], data['interaction'])])

    bin_interaction_cid = {} ####不是零的互作
    for i in bin_cid_interaction.keys():
        if float(bin_cid_interaction[i]) != 0:
            bin_interaction_cid[i] = float(bin_cid_interaction[i])


    bin_name_cid = []  # bin序号的列表
    for e1 in bin_cid_location.keys():
        bin_name_cid.append(e1)

    interaction_frequence_cid = []  # 统计互作频率的列表 # 获取最大的互作频率
    for e1 in bin_interaction_cid.values():
        interaction_frequence_cid.append(float(e1))
    interaction_frequence_cid.sort()
    bin_interaction_iced_cid = {}  # 把之前的互作频率改成0-1的互作频率
    bin_interaction_iced_1_cid = {}
    for e2 in bin_interaction_cid.keys():
        e3 = float(bin_interaction_cid[e2]) / interaction_frequence_cid[-1]  # 这个是将互作频率除最大的互作频率，这样就可以换算成0-1了
        e6 = e2.split(',', 1)
        e7 = e6[0] + ',' + e6[1]
        e8 = e6[1] + ',' + e6[0]
        bin_interaction_iced_cid[e7] = e3
        bin_interaction_iced_cid[e8] = e3  # 因为要涉及到正反的基因对，所以正反我都加上了
        bin_interaction_iced_1_cid[e2] = e3 # 这全是不为零的

    ####3
    # 先把每个bin对应的左边的有信息的十个bin，右边对应的十个bin的列表的字典得到
    # 每个bin对应一个hic scorces就是互作频率，还对应一个得到一个配对t检验的t值
    bin_light_right_cid = {}  # bin和左右对应的bin序号的列表列表的字典
    for r1 in bin_name_cid:  # 先把左右可能想连接的bin的左右bin序号得出来
        light = []
        right = []
        light_right = [] # 先左后右 0：左，1：右
        r2 = 1 # 统计左右互作序号的数量
        while len(light) < 10: # 左右序号没有10个时，就一直进行循环
            r3 = int(r1) - r2 # bin的r1的左边第r2个
            r4 = int(r1) + r2 # bin的r2的右边第r2个
            if r3 < bin_begin_number: #当左边的bin序号小于初始bin序号时，这时就需要变成末尾的bin序号了
                r5 = int(bin_name_cid[-1]) - (bin_begin_number - r3 - 1) # 左边的bin序号
                r6 = r1 + ',' + str(r5)  # 这个就可以得到基因和左边的基因对了，用来判断这个基因对有没有信息
                r7 = int(r1) + r2  # 右边的bin序号
                r8 = r1 + ',' + str(r7)  # 这个就得到的基因右边的基因对
                light.append(r6)
                right.append(r8)
                r2 = r2 + 1
            elif r4 > int(bin_name_cid[-1]): # r4-int(bin_name_cid[-1])是末位置的右边第几个，直接用bin_name_cid的索引进行获取
                il5 = bin_name_cid[(r4 - int(bin_name_cid[-1])) - 1]
                il7 = int(r1) - r2 #左边的序号
                w1 = r1 + ',' + str(il7)
                w2 = r1 + ',' + str(il5)
                light.append(w1)
                right.append(w2)
                r2 = r2 + 1
            else: #既不是在旁边是在里面
                t1 = int(r1) + r2 # 右边的序号
                t2 = int(r1) - r2 # 左边的序号
                t3 = r1 + ',' + str(t2) # 左边的互作bin
                t4 = r1 + ',' + str(t1) # 右边的互作bin
                light.append(t3)
                right.append(t4)
                r2 = r2 + 1
        light_right.append(light)
        light_right.append(right)
        bin_light_right_cid[r1] = light_right

    ####5
    # 得到每个bin和对应的左边的互作bin对应的互作频率的log2列表，和右边的互作bin对应的互作频率的log2列表的字典
    bin_light_right_log_cid = {}
    bin_light_right_interaction_cid = {}  # 得到每个bin和对应的左边的互作bin对应的互作频率的列表，和右边的互作bin对应的互作频率的列表的字典
    for o1, o2 in bin_light_right_cid.items():
        light_interaction = []
        right_interaction = []
        light_right_interaction = []
        light_log = []
        right_log = []
        light_right_log = []
        size_list = 0  # 检测这个bin符合条件的左右交互有几个，如果size_list为0，则不添加这个
        for o3 in range(len(o2[0])):
            if o2[0][o3] in bin_interaction_iced_cid and o2[1][o3] in bin_interaction_iced_cid:  # 两边相同位置同时有值
                if bin_interaction_iced_cid[o2[0][o3]] != 0 and bin_interaction_iced_cid[o2[1][o3]] != 0:  # 两边相同位置同时都不等于0
                    o4 = np.log2(bin_interaction_iced_cid[o2[0][o3]])
                    o5 = np.log2(bin_interaction_iced_cid[o2[1][o3]])
                    light_interaction.append(bin_interaction_iced_cid[o2[0][o3]])
                    right_interaction.append(bin_interaction_iced_cid[o2[1][o3]])
                    light_log.append(o4)
                    right_log.append(o5)
                    size_list = size_list + 1
        if size_list > 4:  # 如果一个符合条件的都没有，就不添加# 最少需要三组，要不然数据不准
            light_right_log.append(light_log)
            light_right_log.append(right_log)
            light_right_interaction.append(light_interaction)
            light_right_interaction.append(right_interaction)
            bin_light_right_log_cid[o1] = light_right_log
            bin_light_right_interaction_cid[o1] = light_right_interaction

    bin_ttest_p_cid = {}  # bin和对应的t值和p值的列表
    for i1, i2 in bin_light_right_log_cid.items():
        i3 = []
        light_l = np.array(i2[0], dtype=float)
        right_r = np.array(i2[1], dtype=float)
        t, p_value = stats.ttest_rel(light_l, right_r)
        i3.append(t)
        i3.append(p_value)
        bin_ttest_p_cid[i1] = i3

    # 先将显著的bin提取出来
    cid_scale = []
    for i in bin_ttest_p_cid.keys():
        if bin_ttest_p_cid[i][0] <= -1.81:
            cid_scale.append(-int(i))
        elif bin_ttest_p_cid[i][0] >= 1.81:
            cid_scale.append(int(i))

    cid_location_number = []  # cid的序号
    j = 0
    for i in cid_scale:
        if (i) * (-1) ** j < 0:
            cid_location_number.append(i)
            j = j + 1

    cid_location_negative = []
    k = bin_begin_number
    for i in range(len(cid_location_number)):
        scale1 = []
        if cid_location_number[i] < 0:
            scale1.append(k)
            scale1.append(-(cid_location_number[i]) - 1)
            k = -cid_location_number[i]
            cid_location_negative.append(scale1)

    cid_location = {}
    for i in range(len(cid_location_negative)):
        scale3 = []
        if i == len(cid_location_negative) - 1:
            scale2 = []
            scale2.append(bin_cid_location[str(cid_location_negative[i][1])][1])
            scale2.append(bin_cid_location[bin_name_cid[-1]][1])
            cid_location[0] = scale2
            scale3.append(bin_cid_location[str(cid_location_negative[i][0])][0])
            scale3.append(bin_cid_location[str(cid_location_negative[i][1])][0])
            cid_location[i + 1] = scale3
        else:
            scale3.append(bin_cid_location[str(cid_location_negative[i][0])][0])
            scale3.append(bin_cid_location[str(cid_location_negative[i][1])][0])
            cid_location[i + 1] = scale3

    return cid_location,bin_ttest_p_cid,cid_location_negative,bin_cid_location

def bin_location_cid(genome_length,resolution_fp,bin_begin_number):
    if (genome_length % resolution_fp) == 0: # 如果没有余值，则不需要除以的商加1
        bin_count_bed = (genome_length // resolution_fp)
    else:
        bin_count_bed = (genome_length // resolution_fp) + 1
    bin_name_bed = []
    for bin_count_bed_list in range(bin_begin_number, bin_count_bed + bin_begin_number):  # bin_number_begin的初始值是0
        bin_name_bed.append(bin_count_bed_list)

    bin_cid_location = {}
    bin_bed_begin = 0
    bin_bed_finna = resolution_fp
    for bin_name_bed_list in bin_name_bed:
        bin_cid_location[str(bin_name_bed_list)] = [bin_bed_begin, bin_bed_finna]
        bin_bed_begin = bin_bed_begin + resolution_fp
        bin_bed_finna = bin_bed_finna + resolution_fp


    bin_name_cid = []  # bin序号的列表
    for e1 in bin_cid_location.keys():
        bin_name_cid.append(e1)

    ####3
    # 先把每个bin对应的左边的有信息的十个bin，右边对应的十个bin的列表的字典得到
    # 每个bin对应一个hic scorces就是互作频率，还对应一个得到一个配对t检验的t值
    bin_light_right_cid = {}  # bin和左右对应的bin序号的列表列表的字典
    for r1 in bin_name_cid:  # 先把左右可能想连接的bin的左右bin序号得出来
        light = []
        right = []
        light_right = []  # 先左后右 0：左，1：右
        r2 = 1  # 统计左右互作序号的数量
        while len(light) < 10:  # 左右序号没有10个时，就一直进行循环
            r3 = int(r1) - r2  # bin的r1的左边第r2个
            r4 = int(r1) + r2  # bin的r2的右边第r2个
            if r3 < bin_begin_number:  # 当左边的bin序号小于初始bin序号时，这时就需要变成末尾的bin序号了
                r5 = int(bin_name_cid[-1]) - (bin_begin_number - r3 - 1)  # 左边的bin序号
                r6 = r1 + ',' + str(r5)  # 这个就可以得到基因和左边的基因对了，用来判断这个基因对有没有信息
                r7 = int(r1) + r2  # 右边的bin序号
                r8 = r1 + ',' + str(r7)  # 这个就得到的基因右边的基因对
                light.append(r6)
                right.append(r8)
                r2 = r2 + 1
            elif r4 > int(bin_name_cid[-1]):  # r4-int(bin_name_cid[-1])是末位置的右边第几个，直接用bin_name_cid的索引进行获取
                il5 = bin_name_cid[(r4 - int(bin_name_cid[-1])) - 1]
                il7 = int(r1) - r2  # 左边的序号
                w1 = r1 + ',' + str(il7)
                w2 = r1 + ',' + str(il5)
                light.append(w1)
                right.append(w2)
                r2 = r2 + 1
            else:  # 既不是在旁边是在里面
                t1 = int(r1) + r2  # 右边的序号
                t2 = int(r1) - r2  # 左边的序号
                t3 = r1 + ',' + str(t2)  # 左边的互作bin
                t4 = r1 + ',' + str(t1)  # 右边的互作bin
                light.append(t3)
                right.append(t4)
                r2 = r2 + 1
        light_right.append(light)
        light_right.append(right)
        bin_light_right_cid[r1] = light_right

    return bin_cid_location,bin_light_right_cid,bin_name_bed[-1],bin_name_bed[0]

def gene_cid_divide(bin_cid_location ,interaction_matrix_filename,bin_light_right_cid,now_path_np,now_seq_np,input_folder,bin_begin_number):
    """
    将bin之间的互作转换为基因间的互作
    :param gene_location_np1: 基因的位置信息的字典
    :param gene_name_np: 基因的名字的形参
    :param gene_gene_interaction_name_np: 互作基因名的列表形参
    :param interaction_matrix_filename: 标准化的bin之间的互作频率文件
    :param resolution_fp: 识别cid的分辨率
    bin_begin_number: 初始bin的序号
    :return: gene_interaction_sort: 基因间的互作频率，从小到大排列；     gene_length: 基因和其对应的长度；    gene1_bin: 基因和对应的bin序号；     gene_interaction: 基因对和对应的互作频率
    """
    import os
    import pandas as pd
    import numpy as np
    from scipy import stats
    ####导入备注文件，得到bin和它的位置信息的字典
    os.chdir(now_path_np + now_seq_np + input_folder)

    ####导入互作数据，得到互作bin序号和互作频率的字典
    data = pd.read_csv(interaction_matrix_filename, sep='\t', header=None)
    bin_cid_interaction = {}
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
            bin_cid_interaction[bin_bin_name_list] = data_matrix.iloc[
                int(bin_bin_name_list1[0]), int(bin_bin_name_list1[1])]
    else:
        data = pd.read_csv(interaction_matrix_filename, sep='\t', header=None)
        sub = ['bin1', 'bin2', 'interaction']  # 赋予列的名字
        data.columns = sub
        # df = data.set_index(['bin1','bin2'])

        # 转为字典
        bin_cid_interaction = dict(
            [(str(i) + ',' + str(a), str(b)) for i, a, b in zip(data['bin1'], data['bin2'], data['interaction'])])

    bin_interaction_cid = {} ####不是零的互作
    for i in bin_cid_interaction.keys():
        if float(bin_cid_interaction[i]) != 0:
            bin_interaction_cid[i] = float(bin_cid_interaction[i])


    bin_name_cid = []  # bin序号的列表
    for e1 in bin_cid_location.keys():
        bin_name_cid.append(e1)

    interaction_frequence_cid = []  # 统计互作频率的列表 # 获取最大的互作频率
    for e1 in bin_interaction_cid.values():
        interaction_frequence_cid.append(float(e1))
    interaction_frequence_cid.sort()
    bin_interaction_iced_cid = {}  # 把之前的互作频率改成0-1的互作频率
    bin_interaction_iced_1_cid = {}
    for e2 in bin_interaction_cid.keys():
        e3 = float(bin_interaction_cid[e2]) / interaction_frequence_cid[-1]  # 这个是将互作频率除最大的互作频率，这样就可以换算成0-1了
        e6 = e2.split(',', 1)
        e7 = e6[0] + ',' + e6[1]
        e8 = e6[1] + ',' + e6[0]
        bin_interaction_iced_cid[e7] = e3
        bin_interaction_iced_cid[e8] = e3  # 因为要涉及到正反的基因对，所以正反我都加上了
        bin_interaction_iced_1_cid[e2] = e3 # 这全是不为零的


    ####5
    # 得到每个bin和对应的左边的互作bin对应的互作频率的log2列表，和右边的互作bin对应的互作频率的log2列表的字典
    bin_light_right_log_cid = {}
    bin_light_right_interaction_cid = {}  # 得到每个bin和对应的左边的互作bin对应的互作频率的列表，和右边的互作bin对应的互作频率的列表的字典
    for o1, o2 in bin_light_right_cid.items():
        light_interaction = []
        right_interaction = []
        light_right_interaction = []
        light_log = []
        right_log = []
        light_right_log = []
        size_list = 0  # 检测这个bin符合条件的左右交互有几个，如果size_list为0，则不添加这个
        for o3 in range(len(o2[0])):
            if o2[0][o3] in bin_interaction_iced_cid and o2[1][o3] in bin_interaction_iced_cid:  # 两边相同位置同时有值
                if bin_interaction_iced_cid[o2[0][o3]] != 0 and bin_interaction_iced_cid[o2[1][o3]] != 0:  # 两边相同位置同时都不等于0
                    o4 = np.log2(bin_interaction_iced_cid[o2[0][o3]])
                    o5 = np.log2(bin_interaction_iced_cid[o2[1][o3]])
                    light_interaction.append(bin_interaction_iced_cid[o2[0][o3]])
                    right_interaction.append(bin_interaction_iced_cid[o2[1][o3]])
                    light_log.append(o4)
                    right_log.append(o5)
                    size_list = size_list + 1
        if size_list > 4:  # 如果一个符合条件的都没有，就不添加# 最少需要三组，要不然数据不准
            light_right_log.append(light_log)
            light_right_log.append(right_log)
            light_right_interaction.append(light_interaction)
            light_right_interaction.append(right_interaction)
            bin_light_right_log_cid[o1] = light_right_log
            bin_light_right_interaction_cid[o1] = light_right_interaction

    bin_ttest_p_cid = {}  # bin和对应的t值和p值的列表
    for i1, i2 in bin_light_right_log_cid.items():
        i3 = []
        light_l = np.array(i2[0], dtype=float)
        right_r = np.array(i2[1], dtype=float)
        t, p_value = stats.ttest_rel(light_l, right_r)
        i3.append(t)
        i3.append(p_value)
        bin_ttest_p_cid[i1] = i3

    # 先将显著的bin提取出来
    cid_scale = []
    for i in bin_ttest_p_cid.keys():
        if bin_ttest_p_cid[i][0] <= -1.81:
            cid_scale.append(-int(i))
        elif bin_ttest_p_cid[i][0] >= 1.81:
            cid_scale.append(int(i))

    cid_location_number = []  # cid的序号
    j = 0
    for i in cid_scale:
        if (i) * (-1) ** j < 0:
            cid_location_number.append(i)
            j = j + 1

    cid_location_negative = []
    k = bin_begin_number
    for i in range(len(cid_location_number)):
        scale1 = []
        if cid_location_number[i] < 0:
            scale1.append(k)
            scale1.append(-(cid_location_number[i]) - 1)
            k = -cid_location_number[i]
            cid_location_negative.append(scale1)

    cid_location = {}
    for i in range(len(cid_location_negative)):
        scale3 = []
        if i == len(cid_location_negative) - 1:
            scale2 = []
            scale2.append(bin_cid_location[str(cid_location_negative[i][1])][1])
            scale2.append(bin_cid_location[bin_name_cid[-1]][1])
            cid_location[0] = scale2
            scale3.append(bin_cid_location[str(cid_location_negative[i][0])][0])
            scale3.append(bin_cid_location[str(cid_location_negative[i][1])][0])
            cid_location[i + 1] = scale3
        else:
            scale3.append(bin_cid_location[str(cid_location_negative[i][0])][0])
            scale3.append(bin_cid_location[str(cid_location_negative[i][1])][0])
            cid_location[i + 1] = scale3

    return cid_location,bin_ttest_p_cid,cid_location_negative

def gene_cid_divide2(bin_cid_location ,interaction_matrix_filename,bin_light_right_cid,now_path_np,now_seq_np,input_folder,bin_begin_number,bin_count_all,resolution_np):
    """
    将bin之间的互作转换为基因间的互作
    :param gene_location_np1: 基因的位置信息的字典
    :param gene_name_np: 基因的名字的形参
    :param gene_gene_interaction_name_np: 互作基因名的列表形参
    :param interaction_matrix_filename: 标准化的bin之间的互作频率文件
    :param resolution_fp: 识别cid的分辨率
    bin_begin_number: 初始bin的序号
    :return: gene_interaction_sort: 基因间的互作频率，从小到大排列；     gene_length: 基因和其对应的长度；    gene1_bin: 基因和对应的bin序号；     gene_interaction: 基因对和对应的互作频率
    """
    import os
    import pandas as pd
    import numpy as np
    from scipy import stats
    ####导入备注文件，得到bin和它的位置信息的字典
    os.chdir(now_path_np + now_seq_np + input_folder)

    ####导入互作数据，得到互作bin序号和互作频率的字典
    data = pd.read_csv(interaction_matrix_filename, sep='\t', header=None)
    bin_cid_interaction = {}
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
            bin_cid_interaction[bin_bin_name_list] = data_matrix.iloc[
                int(bin_bin_name_list1[0]), int(bin_bin_name_list1[1])]
    else:
        data = pd.read_csv(interaction_matrix_filename, sep='\t', header=None)
        sub = ['bin1', 'bin2', 'interaction']  # 赋予列的名字
        data.columns = sub
        # df = data.set_index(['bin1','bin2'])

        # 转为字典
        bin_cid_interaction = dict(
            [(str(i) + ',' + str(a), str(b)) for i, a, b in zip(data['bin1'], data['bin2'], data['interaction'])])

    bin_interaction_cid = {} ####不是零的互作
    for i in bin_cid_interaction.keys():
        if float(bin_cid_interaction[i]) != 0:
            bin_interaction_cid[i] = float(bin_cid_interaction[i])


    bin_name_cid = []  # bin序号的列表
    for e1 in bin_cid_location.keys():
        bin_name_cid.append(e1)

    interaction_frequence_cid = []  # 统计互作频率的列表 # 获取最大的互作频率
    for e1 in bin_interaction_cid.values():
        interaction_frequence_cid.append(float(e1))
    interaction_frequence_cid.sort()
    bin_interaction_iced_cid = {}  # 把之前的互作频率改成0-1的互作频率
    bin_interaction_iced_1_cid = {}
    for e2 in bin_interaction_cid.keys():
        e3 = float(bin_interaction_cid[e2]) / interaction_frequence_cid[-1]  # 这个是将互作频率除最大的互作频率，这样就可以换算成0-1了
        e6 = e2.split(',', 1)
        e7 = e6[0] + ',' + e6[1]
        e8 = e6[1] + ',' + e6[0]
        bin_interaction_iced_cid[e7] = e3
        bin_interaction_iced_cid[e8] = e3  # 因为要涉及到正反的基因对，所以正反我都加上了
        bin_interaction_iced_1_cid[e2] = e3 # 这全是不为零的


    ####5
    # 得到每个bin和对应的左边的互作bin对应的互作频率的log2列表，和右边的互作bin对应的互作频率的log2列表的字典
    bin_light_right_log_cid = {}
    bin_light_right_interaction_cid = {}  # 得到每个bin和对应的左边的互作bin对应的互作频率的列表，和右边的互作bin对应的互作频率的列表的字典
    for o1, o2 in bin_light_right_cid.items():
        light_interaction = []
        right_interaction = []
        light_right_interaction = []
        light_log = []
        right_log = []
        light_right_log = []
        size_list = 0  # 检测这个bin符合条件的左右交互有几个，如果size_list为0，则不添加这个
        for o3 in range(len(o2[0])):
            if o2[0][o3] in bin_interaction_iced_cid and o2[1][o3] in bin_interaction_iced_cid:  # 两边相同位置同时有值
                if bin_interaction_iced_cid[o2[0][o3]] != 0 and bin_interaction_iced_cid[o2[1][o3]] != 0:  # 两边相同位置同时都不等于0
                    o4 = np.log2(bin_interaction_iced_cid[o2[0][o3]])
                    o5 = np.log2(bin_interaction_iced_cid[o2[1][o3]])
                    light_interaction.append(bin_interaction_iced_cid[o2[0][o3]])
                    right_interaction.append(bin_interaction_iced_cid[o2[1][o3]])
                    light_log.append(o4)
                    right_log.append(o5)
                    size_list = size_list + 1
        if size_list > 4:  # 如果一个符合条件的都没有，就不添加# 最少需要三组，要不然数据不准
            light_right_log.append(light_log)
            light_right_log.append(right_log)
            light_right_interaction.append(light_interaction)
            light_right_interaction.append(right_interaction)
            bin_light_right_log_cid[o1] = light_right_log
            bin_light_right_interaction_cid[o1] = light_right_interaction

    bin_ttest_p_cid = {}  # bin和对应的t值和p值的列表
    for i1, i2 in bin_light_right_log_cid.items():
        i3 = []
        light_l = np.array(i2[0], dtype=float)
        right_r = np.array(i2[1], dtype=float)
        t, p_value = stats.ttest_rel(light_l, right_r)
        i3.append(t)
        i3.append(p_value)
        bin_ttest_p_cid[i1] = i3

    # 先将显著的bin提取出来
    cid_scale = []
    for i in bin_ttest_p_cid.keys():
        if bin_ttest_p_cid[i][0] <= -1.81:
            cid_scale.append(-int(i))
        elif bin_ttest_p_cid[i][0] >= 1.81:
            cid_scale.append(int(i))

    cid_location_number = []  # cid的序号
    j = 0
    for i in cid_scale:
        if (i) * (-1) ** j < 0:
            cid_location_number.append(i)
            j = j + 1

    cid_location_negative = []
    k = bin_begin_number
    for i in range(len(cid_location_number)):
        scale1 = []
        if cid_location_number[i] < 0:
            scale1.append(k)
            scale1.append(-(cid_location_number[i]) - 1)
            k = -cid_location_number[i]
            cid_location_negative.append(scale1)

    cid_location = {}
    for i in range(len(cid_location_negative)):
        scale3 = []
        if i == len(cid_location_negative) - 1:
            scale2 = []
            scale2.append(bin_cid_location[str(cid_location_negative[i][1])][1]-(bin_count_all*resolution_np))
            scale2.append(bin_cid_location[bin_name_cid[-1]][1]-(bin_count_all*resolution_np))
            cid_location[0] = scale2
            scale3.append(bin_cid_location[str(cid_location_negative[i][0])][0]-(bin_count_all*resolution_np))
            scale3.append(bin_cid_location[str(cid_location_negative[i][1])][0]-(bin_count_all*resolution_np))
            cid_location[i + 1] = scale3
        else:
            scale3.append(bin_cid_location[str(cid_location_negative[i][0])][0]-(bin_count_all*resolution_np))
            scale3.append(bin_cid_location[str(cid_location_negative[i][1])][0]-(bin_count_all*resolution_np))
            cid_location[i + 1] = scale3

    return cid_location,bin_ttest_p_cid,cid_location_negative

def bin_location_cid2(genome_length,resolution_fp,bin_count_np):
    if (genome_length % resolution_fp) == 0:  # 如果没有余值，则不需要除以的商加1
        bin_count_bed = (genome_length // resolution_fp)
    else:
        bin_count_bed = (genome_length // resolution_fp) + 1
    bin_name_bed = []
    for bin_count_bed_list in range(bin_count_np + 1, bin_count_bed + bin_count_np + 1 ):  # bin_number_begin的初始值是0
        bin_name_bed.append(bin_count_bed_list)

    bin_cid_location = {}
    bin_bed_begin = bin_count_np * resolution_fp
    bin_bed_finna = resolution_fp * (bin_count_np + 1)
    for bin_name_bed_list in bin_name_bed:
        bin_cid_location[str(bin_name_bed_list)] = [bin_bed_begin, bin_bed_finna]
        bin_bed_begin = bin_bed_begin + resolution_fp
        bin_bed_finna = bin_bed_finna + resolution_fp


    bin_name_cid = []  # bin序号的列表
    for e1 in bin_cid_location.keys():
        bin_name_cid.append(e1)

    ####3
    # 先把每个bin对应的左边的有信息的十个bin，右边对应的十个bin的列表的字典得到
    # 每个bin对应一个hic scorces就是互作频率，还对应一个得到一个配对t检验的t值
    bin_light_right_cid = {}  # bin和左右对应的bin序号的列表列表的字典
    for r1 in bin_name_cid:  # 先把左右可能想连接的bin的左右bin序号得出来
        light = []
        right = []
        light_right = [] # 先左后右 0：左，1：右
        r2 = 1 # 统计左右互作序号的数量
        while len(light) < 10:
            r3 = int(r1) - r2  # 这是bin的序号在左边的第几个，主要是判断是不是首尾连接了 #主要是因为原核的染色体是一个圈
            r4 = int(r1) + r2  # 这是bin的序号在左边的第几个，主要是判断是不是首尾连接了
            if r3 < int(bin_name_cid[0]):  # 如果左边的10个小于0的话
                r5 = r3 + int(bin_name_cid[-1]) + 1 - (bin_count_np + 1)  # 左边的bin序号
                r6 = r1 + ',' + str(r5)  # 这个就可以得到基因和左边的基因对了，用来判断这个基因对有没有信息
                r7 = int(r1) + r2  # 右边的bin序号
                r8 = r1 + ',' + str(r7)  # 这个就得到的基因右边的基因对
                light.append(r6)
                right.append(r8)
                r2 = r2 + 1
            elif r4 > int(bin_name_cid[-1]):  # 如果右边的序号大于463时
                t1 = bin_name_cid[(r4 - int(bin_name_cid[-1])) - 1] # 右边的bin序号
                t2 = r1 + ',' + str(t1)  # bin和右边bin的bin对
                t3 = int(r1) - r2  # 左边的bin序号
                t4 = r1 + ',' + str(t3)  # bin和左边bin的bin对
                light.append(t4)
                right.append(t2)
                r2 = r2 + 1
            else:
                y1 = int(r1) + r2  # 右边的bin序号
                y2 = r1 + ',' + str(y1)  # bin和右边的bin序号的bin对
                y3 = int(r1) - r2  # 左边的Bin序号
                y4 = r1 + ',' + str(y3)  # bin和左边的bin序号的bin对
                light.append(y4)
                right.append(y2)
                r2 = r2 + 1
        light_right.append(light)
        light_right.append(right)
        bin_light_right_cid[r1] = light_right

    return bin_cid_location,bin_light_right_cid,bin_name_bed[-1],bin_name_bed[0]