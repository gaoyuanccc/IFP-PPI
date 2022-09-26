def gene_gene_distance(genome_length_fp,gene_location_np):
    """
通过基因的位置信息和基因的长度信息，获取基因对之间的线性距离的字典,并且将基因对按照线性距离从大到小排序
    :param gene_location_np: 基因的位置信息的字典
    :param genome_length_fp: 原核基因组的长度
    :return: gene_distance_sort: 返回的是基因对和对应的线性距离从大到小排序的字典
    gene_name: 基因名的列表
    gene_gene_interaction_name: 互作的基因名，不包括自己和自己的互作的基因
    gene_distance: 互作的基因对名字和对应的线性距离
    """
    ####得到一个基因对和它们之间的线性距离的字典 基因对：线性距离
    # 先得到一个gene名的列表
    gene_name = []
    for l2_list2 in gene_location_np.keys():
        gene_name.append(l2_list2)
    # 再得到一个互作基因名的列表(不包括自己和自己互作的基因对)
    gene_gene_interaction_name = []
    for i_gene_name in range(len(gene_name)):
        i_gene_name_2 = i_gene_name + 1  # 保证比i_gene_name大一个，防止有一样的互作名
        while i_gene_name_2 < len(gene_name):
            i_gene_name_3 = '%s,%s' % (gene_name[i_gene_name], gene_name[i_gene_name_2])
            gene_gene_interaction_name.append(i_gene_name_3)
            i_gene_name_2 = i_gene_name_2 + 1
    # 得比较基因的位置，必须是位置大的末尾位置减去位置小的初始位置 基因长度：4639675
    gene_distance = {}
    for r1 in gene_gene_interaction_name:  # 首先得到基因对
        r2 = r1.split(',', 1)
        if float(gene_location_np[r2[0]][1]) < float(gene_location_np[r2[1]][0]):  # 如果gene1在gene2的前面的话
            r3 = float(gene_location_np[r2[1]][0]) - float(gene_location_np[r2[0]][1])  # gene1和gene2之间的距离
            if r3 <= genome_length_fp // 2:  # 如果基因之间的距离小于一半
                gene_distance[r1] = r3  # 那么基因间的距离就是r3
            else:  # 如果基因的距离大于一半
                r4 = genome_length_fp - float(gene_location_np[r2[1]][1])  # 这是gene2的末尾位置到基因末尾的长度
                r5 = r4 + float(gene_location_np[r2[0]][0])  # 基因间的长度就等于gene1的初始位置加上gene2的末尾位置到基因末尾的长度
                gene_distance[r1] = r5
        elif float(gene_location_np[r2[1]][1]) < float(gene_location_np[r2[0]][0]):  # 如果gene1在gene2的后面的话
            r6 = float(gene_location_np[r2[0]][0]) - float(gene_location_np[r2[1]][1])  # gene2和gene1之间的距离
            if r6 <= genome_length_fp // 2:  # 如果基因之间的距离小于一半
                gene_distance[r1] = r6  # 那么基因间的距离就是r6
            else:  # 如果基因间的距离大于一半
                r7 = genome_length_fp - float(gene_location_np[r2[0]][1])
                r8 = r7 + float(gene_location_np[r2[1]][0])
                gene_distance[r1] = r8
        else:  # 除了上面两种，其他的就是gene1和gene2有交集，所以距离全都归为0
            gene_distance[r1] = 0
    # 把基因对按照线性距离的大小从大到小进行排序,保存为列表

    gene_distance_sort = sorted(gene_distance.items(), key=lambda x: x[1], reverse=True)  # 已经是列表了
    return  gene_name, gene_gene_interaction_name, gene_distance