
def cid_gene(CID_location,gene_location,gene_name,resolution_fp):
    """

    :param CID_location: cid的位置信息字典
    :param gene_location: 基因的位置信息字典
    :param gene_name: 基因名列表
    :param resolution_fp: 识别cids时的分辨率
    :return:
    """
    CID_gene_list = []
    t7 = ['boundary']
    for t1 in CID_location.keys():
        t2 = []
        t2.append(t1)
        for t3 in gene_name:
            if float(CID_location[t1][0]) <= float(gene_location[t3][0]) <= float(CID_location[t1][1]):  # 有一部分在CID中
                t2.append(t3)
            elif CID_location[t1][1] < int(gene_location[t3][0]) < (CID_location[t1][1] + resolution_fp):  # 正好在边界内的
                t7.append(t3)
        CID_gene_list.append(t2)
    CID_gene_list.append(t7)
    CID_gene_1 = []  # 得到一个0和1结合的基因列表
    cid_1_0 = []
    for i in range(len(CID_gene_list)):
        if CID_gene_list[i][0] == 1:
            cid_1_0.append(i)
    for i in range(len(CID_gene_list)):
        if CID_gene_list[i][0] == 0:
            cid_1_0.append(i)
    for i in range(len(CID_gene_list)):
        if CID_gene_list[i][0] != 0 and CID_gene_list[i][0] != 1:
            cid_1_0.append(i)
    h2 = CID_gene_list[cid_1_0[0]] + CID_gene_list[cid_1_0[1]]  # 这是把0和1加起来
    h2.remove(0)  # 然后把列表中的0删除，最后全是CID1的基因了
    CID_gene_1.append(h2)  # 先添加CID1的基因
    for h1 in cid_1_0[2:]:  # 然后从这个列表的第3个开始循环添加
        CID_gene_1.append(CID_gene_list[h1])
    CID_gene = {}  # 是CID和对应的基因的字典
    for g1 in CID_gene_1:
        g3 = []
        for g2 in range(1, len(g1)):
            g3.append(g1[g2])
        CID_gene[str(g1[0])] = g3
    return CID_gene

def cid_gene2(CID_location,gene_location,gene_name,resolution_fp,cid_count_np,chr_number_np):
    """

    :param CID_location: cid的位置信息字典
    :param gene_location: 基因的位置信息字典
    :param gene_name: 基因名列表
    :param resolution_fp: 识别cids时的分辨率
    :return:
    """
    CID_gene_list = []
    t7 = ['boundary']
    for t1 in CID_location.keys():
        t2 = []
        t2.append(t1)
        for t3 in gene_name:
            if float(CID_location[t1][0]) <= float(gene_location[t3][0]) <= float(CID_location[t1][1]):  # 有一部分在CID中
                t2.append(t3)
            elif CID_location[t1][1] < int(gene_location[t3][0]) < (CID_location[t1][1] + resolution_fp):  # 正好在边界内的
                t7.append(t3)
        CID_gene_list.append(t2)
    CID_gene_list.append(t7)
    CID_gene_1 = []  # 得到一个0和1结合的基因列表
    cid_1_0 = []
    for i in range(len(CID_gene_list)):
        if CID_gene_list[i][0] == 1:
            cid_1_0.append(i)
    for i in range(len(CID_gene_list)):
        if CID_gene_list[i][0] == 0:
            cid_1_0.append(i)
    for i in range(len(CID_gene_list)):
        if CID_gene_list[i][0] != 0 and CID_gene_list[i][0] != 1:
            cid_1_0.append(i)
    h2 = CID_gene_list[cid_1_0[0]] + CID_gene_list[cid_1_0[1]]  # 这是把0和1加起来
    h2.remove(0)  # 然后把列表中的0删除，最后全是CID1的基因了
    CID_gene_1.append(h2)  # 先添加CID1的基因
    for h1 in cid_1_0[2:]:  # 然后从这个列表的第3个开始循环添加
        CID_gene_1.append(CID_gene_list[h1])
    CID_gene = {}  # 是CID和对应的基因的字典
    for g1 in CID_gene_1:
        g3 = []
        for g2 in range(1, len(g1)):
            g3.append(g1[g2])
        if g1[0] == 'boundary':
            CID_gene[g1[0] + str(chr_number_np +1)] = g3
        else:
            CID_gene[str(g1[0] + cid_count_np)] = g3
    return CID_gene
