from CGcontent import cal_CGcontent
from PredictStructure import secondary_structure_filtering
import re


def primer_score(list1, conservation_scores):
    primer_seq_score = 0
    nucleotide_count = 0
    nucleotide_total_score = 0

    for k in range(list1[0], list1[1]):
        nucleotide_count += 1
        nucleotide_total_score += conservation_scores[k]

    primer_seq_score = nucleotide_total_score / nucleotide_count

    return primer_seq_score


def get_primers(seq_range, conservation_scores_range, ct, seq, length):
    primers = {}  # primer seq : [primer score, primer start pos, primer end pos]
    for i in range(0, len(seq_range) - 34):

        for j in range(length[0], length[1]+1):

            primer_seq = seq_range[i:i + j]

            if len(re.findall(r'[GC]{6}', primer_seq)) > 0:
                continue

            if '-' in primer_seq:
                continue

            CG_content = cal_CGcontent(primer_seq)

            if 0.3 < CG_content < 0.7:  # CGcontent filtering
                pass
            else:
                continue

            primer_seq_score = primer_score([i, i + j], conservation_scores_range)
            if primer_seq_score < ct:  # conservative filtering
                continue

            primers[primer_seq] = [seq.index(primer_seq), seq.index(primer_seq) + j, primer_seq_score, CG_content]

    primers = secondary_structure_filtering(primers, "DNA")  # secondary structure filtering

    return primers
