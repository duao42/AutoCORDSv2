from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline


def Conserve(infile):
    # Perform multiple sequence alignment using Clustal Omega
    clustalomega_cline = ClustalOmegaCommandline(infile=infile, outfile="alignment.fasta", verbose=True, auto=True)
    stdout, stderr = clustalomega_cline()


def get_seq_conservation_scores(alignment_fasta):
    # Load the aligned sequences
    alignment = AlignIO.read(alignment_fasta, "fasta")
    # Calculate the conservation score for each position in the alignment
    conservation_scores = []
    seq = ""
    for i in range(len(alignment[0])):
        # Slicing a two-dimensional array returns a list of all rows corresponding to the sliced column
        column = alignment[:, i]
        counts = {}
        for nucleotide in column:

            if nucleotide not in counts:
                counts[nucleotide] = 0

            counts[nucleotide] += 1
        keys = [k for k, v in counts.items() if v == max(counts.values())]

        if len(keys) == 1 and keys[0] != "-":
            seq = seq + keys[0]

        else:
            continue
        conservation_score = max(counts.values()) / sum(counts.values())
        conservation_scores.append(conservation_score)

    seq = seq.upper()

    return seq, conservation_scores
