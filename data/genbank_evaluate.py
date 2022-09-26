def genbank_location(genbank_filename,now_path_np,now_seq_np,now_sample_np):
    """
将输入的原核细菌的genbank文件中的基因名和位置信息提取出来
    :param genbank_filename:输入的genbank文件名
    now_sample_np是样品名
    :return: 得到的基因名和对应的位置信息的字典; 以及基因的长度还有基因的方向+1是正的，-1是反的
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    import os
    os.chdir('%s%sinput%s%s%sinteraction_distance' % (now_path_np, now_seq_np, now_seq_np,now_sample_np,now_seq_np))
    # 读取genbank文件
    recs = [rec for rec in SeqIO.parse(genbank_filename, "genbank")]
    for rec in recs: #这一步是读取genbank文件的文件中全是CDS的基因信息
        feats = [feat for feat in rec.features if feat.type == "CDS"]
    gene_location = {}
    gene_direction = {}
    gene_locus_tag = {} #locus_tag名对应的gene名
    join = []
    for i in feats:
        y1 = []
        if len(i.location.parts) == 1:
            y1.append(str(i.location.start).replace('ExactPosition(','').replace(')',''))
            y1.append(str(i.location.end).replace('ExactPosition(','').replace(')',''))
            gene_location[i.qualifiers['locus_tag'][0]] = y1
            gene_direction[i.qualifiers['locus_tag'][0]] = i.location.strand
            if 'gene' in i.qualifiers.keys():
                gene_locus_tag[i.qualifiers['locus_tag'][0]] = i.qualifiers['gene'][0]
        else:
            join.append(i.location)
    genome_length = len(rec.seq)

    return gene_location,genome_length,gene_direction,gene_locus_tag