def kegg_dict(now_path_np,now_seq_np,now_sample_np,kegg_file_name):

    import os
    os.chdir('%s%sinput%s%s%skegg' % (now_path_np, now_seq_np, now_seq_np,now_sample_np,now_seq_np))
    f1 = open(kegg_file_name)
    f2 = f1.readlines()
    f2.pop(0)
    f3 = []
    for i in f2:
        i1 = i.split('\t',-1)
        i2 = []
        i2.append(i1[0])
        i2.append(i1[2])
        f3.append(i2)

    KEGG_dict = {}
    for j in f3:
        j1 = j[1].split(',',-1)
        KEGG_dict[j[0]] = j1

    return KEGG_dict
