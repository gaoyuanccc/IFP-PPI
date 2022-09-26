def string_compare(now_path_np,now_seq_np,string_file_name,now_sample_np):
    """
获取基因对在Stribng数据库中的得分
    :param now_path_np:
    :param now_seq_np:
    :param string_file_name: String数据库文件的名字
    :param gene_pair_np: 筛选到的对应的基因对的列表
    :return: 基因对列表对应的分数
    """
    import os
    os.chdir('%s%sinput%s%s%sstring' % (now_path_np, now_seq_np, now_seq_np,now_sample_np,now_seq_np))
    fstring = open(string_file_name)
    fstring2 = fstring.readlines()
    fstring2.remove(fstring2[0])
    fstring4 = []
    for f2_list in fstring2:
        fstring3 = f2_list.split(' ', -1)
        fstring4.append(fstring3)
    string_ppi_score = []
    for f4_list in fstring4:
        fstring5 = []
        fstring5.append(f4_list[0][-5:].strip())
        fstring5.append(f4_list[1][-5:].strip())
        fstring5.append(f4_list[2].strip())
        fstring5.append(f4_list[3].strip())
        fstring5.append(f4_list[4].strip())
        fstring5.append(f4_list[5].strip())
        fstring5.append(f4_list[6].strip())
        fstring5.append(f4_list[7].strip())
        fstring5.append(f4_list[8].strip())
        fstring5.append(f4_list[9].strip())
        string_ppi_score.append(fstring5)
    string_ppi_score_dict_total = {}
    for string_ppi_score_list in string_ppi_score:
        key_pp = string_ppi_score_list[0] + ',' + string_ppi_score_list[1]
        key_pp1 = string_ppi_score_list[1] + ',' + string_ppi_score_list[0]
        string_ppi_score_dict_total[key_pp] = string_ppi_score_list[9]
        string_ppi_score_dict_total[key_pp1] = string_ppi_score_list[9]
    return string_ppi_score_dict_total

def string_au(now_path_np,now_seq_np,string_file_name,gene_pair_dict,now_sample_np,threshold_np):
    """

    :param now_path_np:
    :param now_seq_np:
    :param string_file_name:
    :param gene_pair_dict: 基因对和对应的交互频率或者线性距离值的字典
    :param now_sample_np:
    """
    import os
    os.chdir('%s%sinput%s%s%sstring' % (now_path_np, now_seq_np, now_seq_np, now_sample_np, now_seq_np))
    fstring = open(string_file_name)
    fstring2 = fstring.readlines()
    fstring2.remove(fstring2[0])
    fstring4 = []
    for f2_list in fstring2:
        fstring3 = f2_list.split(' ', -1)
        fstring4.append(fstring3)
    string_ppi_score = []
    for f4_list in fstring4:
        fstring5 = []
        fstring5.append(f4_list[0][-5:].strip())
        fstring5.append(f4_list[1][-5:].strip())
        fstring5.append(f4_list[2].strip())
        fstring5.append(f4_list[3].strip())
        fstring5.append(f4_list[4].strip())
        fstring5.append(f4_list[5].strip())
        fstring5.append(f4_list[6].strip())
        fstring5.append(f4_list[7].strip())
        fstring5.append(f4_list[8].strip())
        fstring5.append(f4_list[9].strip())
        string_ppi_score.append(fstring5)
    string_ppi_score_dict_total = {}
    for string_ppi_score_list in string_ppi_score:
        key_pp = string_ppi_score_list[0] + ',' + string_ppi_score_list[1]
        key_pp1 = string_ppi_score_list[1] + ',' + string_ppi_score_list[0]
        string_ppi_score_dict_total[key_pp] = string_ppi_score_list[9]
        string_ppi_score_dict_total[key_pp1] = string_ppi_score_list[9]

    gene_true = []
    gene_list = []
    for i1,j1 in gene_pair_dict.items():
        if i1 in string_ppi_score_dict_total.keys():
            if float(string_ppi_score_dict_total[i1]) >= threshold_np:
                gene_true.append(1)
                gene_list.append(j1)
            else:
                gene_true.append(0)
                gene_list.append(j1)
        else:
            gene_true.append(0)
            gene_list.append(j1)

    return gene_true,gene_list

def string_au_distance(now_path_np,now_seq_np,string_file_name,gene_pair_dict,now_sample_np,threshold_np):
    """
    :param now_path_np:
    :param now_seq_np:
    :param string_file_name:
    :param gene_pair_dict: 基因对和对应的交互频率或者线性距离值的字典
    :param now_sample_np:
    """
    import os
    os.chdir('%s%sinput%s%s%sstring' % (now_path_np, now_seq_np, now_seq_np, now_sample_np, now_seq_np))
    fstring = open(string_file_name)
    fstring2 = fstring.readlines()
    fstring2.remove(fstring2[0])
    fstring4 = []
    for f2_list in fstring2:
        fstring3 = f2_list.split(' ', -1)
        fstring4.append(fstring3)
    string_ppi_score = []
    for f4_list in fstring4:
        fstring5 = []
        fstring5.append(f4_list[0][-5:].strip())
        fstring5.append(f4_list[1][-5:].strip())
        fstring5.append(f4_list[2].strip())
        fstring5.append(f4_list[3].strip())
        fstring5.append(f4_list[4].strip())
        fstring5.append(f4_list[5].strip())
        fstring5.append(f4_list[6].strip())
        fstring5.append(f4_list[7].strip())
        fstring5.append(f4_list[8].strip())
        fstring5.append(f4_list[9].strip())
        string_ppi_score.append(fstring5)
    string_ppi_score_dict_total = {}
    for string_ppi_score_list in string_ppi_score:
        key_pp = string_ppi_score_list[0] + ',' + string_ppi_score_list[1]
        key_pp1 = string_ppi_score_list[1] + ',' + string_ppi_score_list[0]
        string_ppi_score_dict_total[key_pp] = string_ppi_score_list[9]
        string_ppi_score_dict_total[key_pp1] = string_ppi_score_list[9]

    gene_true = []
    gene_list = []
    for i1,j1 in gene_pair_dict.items():
        if i1 in string_ppi_score_dict_total.keys():
            if float(string_ppi_score_dict_total[i1]) >= threshold_np:
                gene_true.append(1)
                gene_list.append(2320813-j1)
            else:
                gene_true.append(0)
                gene_list.append(2320813-j1)
        else:
            gene_true.append(0)
            gene_list.append(2320813-j1)
    return gene_true, gene_list

def string_au_list(now_path_np,now_seq_np,string_file_name,gene_pair_list,now_sample_np,threshold_np,gene_pair):
    """

    :param now_path_np:
    :param now_seq_np:
    :param string_file_name:
    :param gene_pair_list:
    :param now_sample_np:
    :param threshold_np:
    :param gene_pair: 基因对和对应的交互频率和线性距离
    """
    import os
    os.chdir('%s%sinput%s%s%sstring' % (now_path_np, now_seq_np, now_seq_np, now_sample_np, now_seq_np))
    fstring = open(string_file_name)
    fstring2 = fstring.readlines()
    fstring2.remove(fstring2[0])
    fstring4 = []
    for f2_list in fstring2:
        fstring3 = f2_list.split(' ', -1)
        fstring4.append(fstring3)
    string_ppi_score = []
    for f4_list in fstring4:
        fstring5 = []
        fstring5.append(f4_list[0][-5:].strip())
        fstring5.append(f4_list[1][-5:].strip())
        fstring5.append(f4_list[2].strip())
        fstring5.append(f4_list[3].strip())
        fstring5.append(f4_list[4].strip())
        fstring5.append(f4_list[5].strip())
        fstring5.append(f4_list[6].strip())
        fstring5.append(f4_list[7].strip())
        fstring5.append(f4_list[8].strip())
        fstring5.append(f4_list[9].strip())
        string_ppi_score.append(fstring5)
    string_ppi_score_dict_total = {}
    for string_ppi_score_list in string_ppi_score:
        key_pp = string_ppi_score_list[0] + ',' + string_ppi_score_list[1]
        key_pp1 = string_ppi_score_list[1] + ',' + string_ppi_score_list[0]
        string_ppi_score_dict_total[key_pp] = string_ppi_score_list[9]
        string_ppi_score_dict_total[key_pp1] = string_ppi_score_list[9]

    gene_true = []
    gene_list = []
    for i1 in gene_pair_list:
        if i1 in string_ppi_score_dict_total.keys():
            if float(string_ppi_score_dict_total[i1]) >= threshold_np:
                gene_true.append(1)
                gene_list.append(gene_pair[i1])
            else:
                gene_true.append(0)
                gene_list.append(gene_pair[i1])
        else:
            gene_true.append(0)
            gene_list.append(gene_pair[i1])
    return gene_true,gene_list

def string_au_list_distance(now_path_np,now_seq_np,string_file_name,gene_pair_list,now_sample_np,threshold_np,gene_pair):
    """

    :param now_path_np:
    :param now_seq_np:
    :param string_file_name:
    :param gene_pair_list:
    :param now_sample_np:
    :param threshold_np:
    :param gene_pair: 基因对和对应的交互频率和线性距离
    """
    import os
    os.chdir('%s%sinput%s%s%sstring' % (now_path_np, now_seq_np, now_seq_np, now_sample_np, now_seq_np))
    fstring = open(string_file_name)
    fstring2 = fstring.readlines()
    fstring2.remove(fstring2[0])
    fstring4 = []
    for f2_list in fstring2:
        fstring3 = f2_list.split(' ', -1)
        fstring4.append(fstring3)
    string_ppi_score = []
    for f4_list in fstring4:
        fstring5 = []
        fstring5.append(f4_list[0][-5:].strip())
        fstring5.append(f4_list[1][-5:].strip())
        fstring5.append(f4_list[2].strip())
        fstring5.append(f4_list[3].strip())
        fstring5.append(f4_list[4].strip())
        fstring5.append(f4_list[5].strip())
        fstring5.append(f4_list[6].strip())
        fstring5.append(f4_list[7].strip())
        fstring5.append(f4_list[8].strip())
        fstring5.append(f4_list[9].strip())
        string_ppi_score.append(fstring5)
    string_ppi_score_dict_total = {}
    for string_ppi_score_list in string_ppi_score:
        key_pp = string_ppi_score_list[0] + ',' + string_ppi_score_list[1]
        key_pp1 = string_ppi_score_list[1] + ',' + string_ppi_score_list[0]
        string_ppi_score_dict_total[key_pp] = string_ppi_score_list[9]
        string_ppi_score_dict_total[key_pp1] = string_ppi_score_list[9]

    gene_true = []
    gene_list = []
    for i1 in gene_pair_list:
        if i1 in string_ppi_score_dict_total.keys():
            if float(string_ppi_score_dict_total[i1]) >= threshold_np:
                gene_true.append(1)
                gene_list.append(2320813 - gene_pair[i1])
            else:
                gene_true.append(0)
                gene_list.append(2320813 - gene_pair[i1])
        else:
            gene_true.append(0)
            gene_list.append(2320813-gene_pair[i1])

    return gene_true,gene_list

def string_experiment_database(now_path_np,now_seq_np,string_file_name,now_sample_np):
    """
只加实验和精准数据库的得分,返回基因对和得分，和分值的列表
    :param now_path_np:
    :param now_seq_np:
    :param string_file_name: String数据库文件的名字
    :param gene_pair_np: 筛选到的对应的基因对的列表
    :return: 基因对列表对应的分数
    """
    import os
    os.chdir('%s%sinput%s%s%sstring' % (now_path_np, now_seq_np, now_seq_np,now_sample_np,now_seq_np))
    f1 = open(string_file_name)
    f2 = f1.readlines()
    f2.pop(0)
    f3 = []
    for i in f2:
        i1 = []
        j = i.split(' ', -1)
        i1.append(j[0].replace('511145.', ''))
        i1.append(j[1].replace('511145.', ''))
        i2 = int(j[6])
        i1.append(i2)
        f3.append(i1)

    string_dict = {}
    for q in f3:
        q1 = q[0] + ',' + q[1]
        q2 = q[1] + ',' + q[0]
        string_dict[q1] = q[2]
        string_dict[q2] = q[2]

    score_list = []  # 得分列表
    for e, r in string_dict.items():
        score_list.append(r)

    return string_dict,score_list
