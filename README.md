# AutoCORDSv2
An automated tool for designing Cas12a crRNA and RPA primers for CORDSv2 (Cas12a-based On-site and Rapid Detection System) Runs on a Linux system.
# Usage
Download the AutoCORDSv2 folder and run the python file named AutoCORDSv2.py.

A simple usage is "python3 -i "

## Dependencies

Python v3.7

biopython v1.79

clustalo v1.2.4

RNAfold v2.5.1

## Argument

options:

  -h, --help            show this help message and exit
  
  -i Input              Path of input file
  
  -t Type               fasta/MSA_fasta.  Multiple genome files or Multiple sequence alignment files generated by Clustal Omega（The default name of MSA file that generated by this programe is ’alignment.fasta‘）
                        
  -cl crRNA Length      int.  The length of crRNA, default = 20.
  
  -pl the Range of Primer Length
                        str.  The range of primer length, default = 28,35.
                        
  -cct crRNA Conserve Threshold
                        float.  default value is 0.99.
                        
  -pct Primer Conserve Threshold
                        float.  default value is 0.99.
                        
  -pst Primer Specificity Threshold
                        float.  The minimum ratio of mismatched bases allowed between the primer and the background genome, default value is 0.3.
                        
  -n N                  The number of processes, default=8.
  
  -range RANGE          The range of conversed seq to search crRNA&primer.
  
  --bgfile BGFILE [BGFILE ...]
                        Path of background fasta file(one or multiple files separated by spaces) for the Primer or crRNA specificity screening.
                        

# Citations

Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3 https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878

Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539 doi:10.1038/msb.2011.75

Camacho, C., Coulouris, G., Avagyan, V. et al. BLAST+: architecture and applications. BMC Bioinformatics 10, 421 (2009). https://doi.org/10.1186/1471-2105-10-421

Lorenz, R., Bernhart, S.H., Höner zu Siederdissen, C. et al. ViennaRNA Package 2.0. Algorithms Mol Biol 6, 26 (2011). https://doi.org/10.1186/1748-7188-6-26
