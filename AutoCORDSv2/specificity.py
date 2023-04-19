def get_hanming(seq1, seq2):
    hanming_dist = 0
    for (c1, c2) in zip(seq1, seq2):
        if c1 != c2:
            hanming_dist += 1
    return hanming_dist


def get_conserved_seq(string, scores, t):
    high_score_chars = []
    current_high_score_chars = ''

    for i in range(len(string)):
        if scores[i] >= t:
            current_high_score_chars += string[i]

        else:
            if len(current_high_score_chars) >= 25:
                high_score_chars.append(current_high_score_chars)
                current_high_score_chars = ''

    if len(current_high_score_chars) >= 25:
        high_score_chars.append(current_high_score_chars)

    return high_score_chars


def str_split_n(string, n):
    split_list = []
    num_splits = len(string) - n +1
    split_list = [string[i:i + n] for i in range(num_splits)]
    return split_list


def reverse_complement(sequence):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    reverse_seq = sequence[::-1]
    reverse_complement_seq = "".join(complement.get(base, base) for base in reverse_seq)
    return reverse_complement_seq