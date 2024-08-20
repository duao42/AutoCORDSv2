import re
from CGcontent import cal_CGcontent
from PredictStructure import secondary_structure_filtering


def find_position(target_str, source_str):
    pattern = re.compile(target_str)
    match_obj = pattern.search(source_str)
    if match_obj:
        return match_obj.start(), match_obj.end()
    else:
        return -1


def crRNA_score(position_list, conservation_scores):
    crRNA_seq_score = 0
    nucleotide_count = 0
    nucleotide_total_score = 0

    for k in range(position_list[0], position_list[1]):
        nucleotide_count += 1
        nucleotide_total_score += conservation_scores[k]

    crRNA_seq_score = nucleotide_total_score / nucleotide_count

    return crRNA_seq_score


def generate_sequences(rule):
    replacements = {"N": ['A', 'T', 'C', 'G'], "R": ['A', 'G'], "Y": ['C', 'T'], "M": ['A', 'C'], "K": ['G', 'T'],
                    "H": ['A', 'C', 'T'], "D": ['A', 'G', 'T'], "B": ['C', 'G', 'T'], "V": ['A', 'C', 'G']}

    def generate_recursive(current, index):
        # 如果当前索引等于规则长度，说明已经生成了一个完整的字符串
        if index == len(rule):
            results.append(current)
            return

        char = rule[index]

        if char in replacements:
            for replacement in replacements[char]:
                generate_recursive(current + replacement, index + 1)
        else:
            generate_recursive(current + char, index + 1)

    results = []
    generate_recursive("", 0)
    return results


def get_crRNAs(seq, PAM, cl, strand):

    pam_list = generate_sequences(PAM)
    crRNAs = []

    # search spacer by pam
    for i in range(0, len(seq)):
        if seq[i:i + len(PAM)] in pam_list and strand == '-':
            if i - cl >= 0 and i + len(PAM) < len(seq):
                crRNAs.append(seq[i - cl: i + len(PAM)])

        if seq[i:i + len(PAM)] in pam_list and strand == '+':
            if i >= 0 and i + len(PAM) + cl < len(seq):
                crRNAs.append(seq[i:i + len(PAM) + cl])

    crRNA_dict = {}
    for crRNA in crRNAs:
        crRNA_dict[crRNA] = []
    return crRNA_dict


def crRNA_1st_filtering_pos(crRNA_dict, seq, conservation_scores, ct):
    # find crRNA position
    crRNA_dict_1 = {}
    for crRNA in crRNA_dict:
        crRNA_dict_1[crRNA] = [(find_position(crRNA, seq)[0]), (find_position(crRNA, seq)[1])]

    # crRNA conserved
    crRNA_dict_2 = {}
    for crRNA in crRNA_dict:
        position_list = crRNA_dict_1[crRNA]
        crRNA_cons_score = crRNA_score(position_list, conservation_scores)
        crRNA_dict_2[crRNA] = crRNA_dict_1[crRNA] + [crRNA_cons_score]

    # crRNA CG content
    crRNA_dict_3 = {}
    for crRNA in crRNA_dict:
        CG_content = cal_CGcontent(crRNA[4:])
        crRNA_dict_3[crRNA] = crRNA_dict_2[crRNA] + [CG_content]

    crRNA_dict_new = {}
    for crRNA in crRNA_dict_3:
        if crRNA_dict_3[crRNA][2] >= ct and 0.3 < crRNA_dict_3[crRNA][3] < 0.7:
            crRNA_dict_new[crRNA] = crRNA_dict_3[crRNA]

    crRNA_dict_new = secondary_structure_filtering(crRNA_dict_new, "RNA_pos")

    return crRNA_dict_new


def crRNA_1st_filtering_rev(crRNA_dict, seq, conservation_scores, ct, cl):
    # find crRNA position
    crRNA_dict_1 = {}
    for crRNA in crRNA_dict:
        crRNA_dict_1[crRNA] = [(find_position(crRNA, seq)[0]), (find_position(crRNA, seq)[1])]

    # crRNA conserved
    crRNA_dict_2 = {}
    for crRNA in crRNA_dict:
        position_list = crRNA_dict_1[crRNA]
        crRNA_cons_score = crRNA_score(position_list, conservation_scores)
        crRNA_dict_2[crRNA] = crRNA_dict_1[crRNA] + [crRNA_cons_score]

    # crRNA CG content
    crRNA_dict_3 = {}
    for crRNA in crRNA_dict:
        CG_content = cal_CGcontent(crRNA[0:cl])
        crRNA_dict_3[crRNA] = crRNA_dict_2[crRNA] + [CG_content]

    crRNA_dict_new = {}
    for crRNA in crRNA_dict_3:
        if crRNA_dict_3[crRNA][2] >= ct and 0.3 < crRNA_dict_3[crRNA][3] < 0.7:
            crRNA_dict_new[crRNA] = crRNA_dict_3[crRNA]

    crRNA_dict_new = secondary_structure_filtering(crRNA_dict_new, "RNA_rev")

    return crRNA_dict_new
