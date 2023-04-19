import subprocess
import os
from specificity import reverse_complement


def secondary_structure_filtering(Primer_or_crRNA_dict, type):
    if type == "DNA":

        with open("RNAfold_primer.fa", "w") as f:

            for key in Primer_or_crRNA_dict:
                f.write(">" + key + '\n' + key + '\n')

        RNAfold_cmd = ["RNAfold -P DNA --noPS --noconv -i RNAfold_primer.fa > RNAfold_primer.txt"]
        process = subprocess.Popen(RNAfold_cmd, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.DEVNULL)  # shell=True ensure the workspace is right
        output = process.communicate()  # wait for subprocess end

        with open("RNAfold_primer.txt", "r") as f:
            n = 1

            for line in f:
                if n % 3 == 2:
                    key = line.strip()

                if n % 3 == 0:
                    min_energy = float(line.strip()[-6:].replace(")", "").replace(" ", ""))
                    if min_energy < -5:
                        del Primer_or_crRNA_dict[key]
                    else:
                        Primer_or_crRNA_dict[key].append(min_energy)

                n += 1

        os.remove('RNAfold_primer.fa')
        os.remove('RNAfold_primer.txt')
        return Primer_or_crRNA_dict

    if type == "RNA_pos":
        with open("RNAfold_crRNA.fa", "w") as f:

            for key in Primer_or_crRNA_dict:
                f.write(">" + key[4:] + '\n' + key + '\n')

        RNAfold_cmd = ["RNAfold --noPS -i RNAfold_crRNA.fa > RNAfold_crRNA.txt"]
        process = subprocess.Popen(RNAfold_cmd, shell=True,
                                   stdout=subprocess.PIPE)  # shell=True ensure the workspace is right
        output = process.communicate()  # wait for subprocess end
        with open("RNAfold_crRNA.txt", "r") as f:
            n = 1
            for line in f:

                if n % 3 == 2:
                    key = line.strip().replace('U', 'T')

                if n % 3 == 0:
                    min_energy = float(line.strip()[-6:].replace(")", "").replace(" ", ""))
                    if min_energy < -5:
                        del Primer_or_crRNA_dict[key]
                    else:
                        Primer_or_crRNA_dict[key].append(min_energy)

                n += 1

        os.remove('RNAfold_crRNA.fa')
        os.remove('RNAfold_crRNA.txt')
        return Primer_or_crRNA_dict

    if type == "RNA_rev":
        with open("RNAfold_crRNA.fa", "w") as f:

            for key in Primer_or_crRNA_dict:
                f.write(">" + reverse_complement(key)[4:] + '\n' + key + '\n')

        RNAfold_cmd = ["RNAfold --noPS -i RNAfold_crRNA.fa > RNAfold_crRNA.txt"]
        process = subprocess.Popen(RNAfold_cmd, shell=True,
                                   stdout=subprocess.PIPE)  # shell=True ensure the workspace is right
        output = process.communicate()  # wait for subprocess end
        with open("RNAfold_crRNA.txt", "r") as f:
            n = 1
            for line in f:

                if n % 3 == 2:
                    key = line.strip().replace('U', 'T')

                if n % 3 == 0:
                    min_energy = float(line.strip()[-6:].replace(")", "").replace(" ", ""))
                    if min_energy < -5:
                        del Primer_or_crRNA_dict[key]
                    else:
                        Primer_or_crRNA_dict[key].append(min_energy)

                n += 1

        os.remove('RNAfold_crRNA.fa')
        os.remove('RNAfold_crRNA.txt')
        return Primer_or_crRNA_dict

