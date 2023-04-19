# Required files

AutoCORDSv2 code

"12L_MSA.fasta" in the Demo folder

background file that download from "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz" and unzip

# Run

Make sure all required files are in the same folder and required dependencies are installed.

Go into that folder and run the following command:

python3 -i 12L_MSA.fasta -t MSA_fasta -cl 20 -pl 30,35 -cct 0.985 -pct 0.985 -pst 0.3 -n 20 --bgfile bg_file.fna

The output can then be seen in the 12L_MSA folder.