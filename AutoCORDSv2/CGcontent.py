def cal_CGcontent(seq):
    CG_count = 0
    seq = seq.upper()

    for str in seq:
        if str == "C" or str == "G":
            CG_count += 1
    return CG_count/len(seq)