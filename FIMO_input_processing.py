import os
from additional_code.read_MANE import *
import random
import numpy as np
from matplotlib import pyplot as plt

dict_MANE = dict_MANE
###############################################################################################################################################
# FUNCTIONS THAT EXTRACT 1. PPMS FROM ATtRACT, 2. PPMS FROM HTSELEX
###############################################################################################################################################


def extract_attract_ppm(attract_ppm):

  dict_attract_ppms = {}
  with open(attract_ppm, "r") as f:
    for line in f:
        line = line.rstrip("\n").split("\t")

        if len(line) == 2:
            counter = 0
            matrix_id = line[0][1:]
            length = int(line[1])
            dict_curr_PPM = {}
            continue

        dict_curr_PPM[counter] = [float(x) for x in line]
        counter += 1
        if counter == length:
            dict_attract_ppms[matrix_id] = dict_curr_PPM

    return dict_attract_ppms

def extract_selex_ppm(selex_ppm, f_only_monomer=True):
    dict_htselex_ppms = {}

    with open(selex_ppm, "r") as f:
        for line in f:
            line = line.rstrip("\n").split("\t")

            if line[0].startswith(">"):
                counter = 0
                gene_id = line[0][1:]
                length = int(line[1])
                mer_status = line[3]
                dict_curr_PPM = {}
                continue

            dict_curr_PPM[counter] = [float(x) for x in line]
            counter += 1
            if counter == length:
                if f_only_monomer and mer_status == "dimeric":
                    continue
                elif gene_id not in dict_MANE:
                    continue
                else:
                    dict_htselex_ppms[gene_id] = dict_curr_PPM


    return dict_htselex_ppms

data_path = os.path.abspath("DATA")

attract_ppms = extract_attract_ppm(os.path.join(data_path, "ATtRACT_ppm.txt"))
htselex_ppms = extract_selex_ppm(os.path.join(data_path, "SELEX2020_ppm.txt"))

#print(">>> DONE EXTRACTING PPMS FROM ATTRACT AND HT-SELEX\n")

###############################################################################################################################################
# FUNCTION GETTING NESTED DICT FROM ATtRACT DATABSE:
# EXPERIMENT:GENE_ID and GENE_ID:SEQUENCE
# + ONLY INCLUDES GENES THAT ARE AVAILABLE IN MANE DATABASE

# IMPORTANT VARIABLE: attract_ppms
###############################################################################################################################################


def read_db_attract(experiment_name, f_organism, f_experiment, f_score=True, f_mutated=False):
    dict_attract_db = {}
    data_path = os.path.abspath("DATA")
    attract_path = os.path.join(os.path.join(data_path, "ATtRACT_db.txt"))

    dict_ppms = extract_attract_ppm(os.path.join(data_path, "ATtRACT_ppm.txt"))
    gene_ids_with_duplicates = []
    with open(attract_path, "r") as dbfile:
        header_attract = dbfile.readline()
        for line in dbfile:
            line = line.split("\t")

            # Define database structure:
            Gene_name = line[0]
            Gene_id = line[1]
            Mutated = line[2]
            Organism = line[3]
            Motif = line[4]
            # Len = line[5]
            Experiment_description = line[6]
            Matrix_id = line[11]
            Score = line[12]

            # Filter for variables defined above:
            if f_organism:
                if Organism != f_organism:
                    continue
            if f_experiment:
                if Experiment_description != f_experiment:
                    continue
            if f_score:
                if "**" not in Score:
                    continue

            # Optional filters:

            if f_mutated:
                if "yes" in Mutated:
                    continue

            if Gene_id not in dict_MANE:
                continue

            if f_experiment == "RNAcompete":
                if Gene_id in gene_ids_with_duplicates:
                    counter_for_matrix_duplicates = gene_ids_with_duplicates.count(Gene_id)
                    gene_ids_with_duplicates.append(Gene_id)
                    Gene_id = Gene_id+"_"+str(counter_for_matrix_duplicates)
                    dict_attract_db[Gene_id] = Matrix_id
                    continue

            gene_ids_with_duplicates.append(Gene_id)

            dict_attract_db[Gene_id] = Matrix_id

        attract_ppms = {}
        for gene_id, matrix_id in dict_attract_db.items():
            attract_ppms[gene_id] = dict_ppms[matrix_id]

    return attract_ppms

attract_ppms = {}
attract_ppms["RNAcomp"] = read_db_attract("RNAcompete", "Homo_sapiens", "RNAcompete", f_score=True, f_mutated=False)
attract_ppms["SELEX"] = read_db_attract("SELEX", "Homo_sapiens", "SELEX", f_score=True, f_mutated=False)

RNAcomp_ppm_num = len(attract_ppms["RNAcomp"].keys())
SELEX_ppm_num = len(attract_ppms["SELEX"].keys())
htselex_ppm_num = len(htselex_ppms.keys())

#print(">>> DONE CREATING PPM DICTIONARIES FOR FURTHER PROCESSING\n")
#print(f">>> NUM OF ATTRACT RNAcomp MOTIFS (pre-filter): {RNAcomp_ppm_num}\n")
#print(f">>> NUM OF ATTRACT RNAcomp MOTIFS (pre-filter): {SELEX_ppm_num}\n")
#print(f">>> NUM OF HT-SELEX MOTIFS (pre-filter): {htselex_ppm_num}\n")

def plot_motif_length_distribution(list_of_motif_lengths):
    list_of_motif_lengths = sorted(list_of_motif_lengths)
    min = np.min(list_of_motif_lengths)
    max = np.max(list_of_motif_lengths)

    fig, ax = plt.subplots( )
    ax.hist(list_of_motif_lengths, bins=np.arange(min, max, 1))

    plt.show( )


#plot_motif_length_distribution(array)





###############################################################################################################################################
# FILTER TRANSCRIPTOME DATABASE FOR UTR3, UTR5, CDS AND FULL TRANSCRIPT
###############################################################################################################################################


subseqs = ["UTR3", "UTR5", "CDS", "cDNA"]
MANE_transcriptome = {}

for id, content in dict_MANE.items():
    MANE_transcriptome[id] = {}
    for subseq in subseqs:
        MANE_transcriptome[id][subseq] = dict_MANE[id][subseq]

#print(">>> DONE SPLITTING MANE SEQUENCES INTO SUBSEQUENCES\n")


def rename_MANE_seqs_duplicated_motifs(MANE_transcriptome, rnacompete_ppms):
    for ppm in rnacompete_ppms:
        if ppm not in MANE_transcriptome:
            ppm_name_without_extension = str(ppm)[:-2]
            MANE_transcriptome[ppm] = MANE_transcriptome[ppm_name_without_extension]

    for key in MANE_transcriptome.keys():
        if "cDNA" in MANE_transcriptome[key].keys():
            MANE_transcriptome[key]["transcript"] = MANE_transcriptome[key].pop('cDNA')
        elif "transcript" in MANE_transcriptome[key]:
            continue

    return MANE_transcriptome

MANE_transcriptome = rename_MANE_seqs_duplicated_motifs(MANE_transcriptome, attract_ppms["RNAcomp"])

###############################################################################################################################################
# FILTER UTRs AND CDS FOR EMPTY/SHORT SEQUENCES;
# DETERMINE A CUTOFF FOR SEQUENCE LENGTH 20
###############################################################################################################################################
def len_cutoff_subseq(cutoff=20):

    for id in MANE_transcriptome.keys():
        for subseq in MANE_transcriptome[id].keys():
            if len(MANE_transcriptome[id][subseq]) <= cutoff:
                MANE_transcriptome[id][subseq] = ""

    #print(">>> DONE FILTERING FOR SUBSEQUENCE LENGTH\n")

len_cutoff_subseq(cutoff=0)




###############################################################################################################################################
# CREATE RANDOMLY SHUFFLED SEQUENCES FROM DINUCLEOTIDE FREQUENCIES
###############################################################################################################################################

def generate_random_di_nuc_seq(length):
    di_data = ({
    "A": [0.29964487076119456, 0.19291583797335501, 0.28081683230189497, 0.22662245896355537],
    "C": [0.30514068401513184, 0.29142157085594184, 0.10437582930948616, 0.29906191581944036],
    "G": [0.2738100837255303, 0.24900349340170375, 0.27287036424103384, 0.2043160586317321],
    "U": [0.1787865306201401, 0.2254019334639803, 0.3002818262111331, 0.29552970970474646]},
    [0.26366181575547193, 0.23828855841092605, 0.24192231774690498, 0.25612730808669704])

    dep_P, starts = di_data
    random_seq = []
    for i in range(length):
        if i == 0:
            random_seq.append(np.random.choice(["A","C","G","U"], size=None, replace=True, p=starts))
        else:
            random_seq.append(np.random.choice(["A","C","G","U"], size=None, replace=True, p=dep_P[random_seq[-1]]))
    return "".join(random_seq)



###############################################################################################################################################
# CREATE RANDOMLY SHUFFLED SEQUENCES IN DIFFERENT WAYS:
# - BY SHUFFLING A MANE TRANSCRIPT
# - BY SHUFFLING THE AUTOLOGOUS TRANSCRIPT
# - BY GENERATING NEW RANDOM SEQUENCES VIA DINUCLEOTIDE FREQUENCIES
###############################################################################################################################################

def get_randomized_sequence_list(transcript_dict, seq, seq_amount, seq_type_key, shuffle_mode="transcriptome"):
    # If transcriptome is used, get sequences from MANE with at least the length of the autologous
    # sequence and save them in list:
    if shuffle_mode == "transcriptome":
        list_usable_seqs = []
        for protein in transcript_dict.keys():
            if len(transcript_dict[protein][seq_type_key]) >= len(seq):
                list_usable_seqs.append(transcript_dict[protein][seq_type_key])

        # Get random sequences from the list until defined amount is met:
        if len(list_usable_seqs) >= seq_amount:
            list_amount = random.sample(list_usable_seqs, seq_amount)
        else:
            list_amount = random.choices(list_usable_seqs, k=seq_amount)
            #print(f"\nNot enough background for a {seq_type_key}, used repetitions!\n")

    # Create randomized background with length of sequence and amount of randomized sequences
    # given by user:
    randomized_sequences = []
    for num in range(seq_amount):
        # Shuffle input sequence to get random sequences:
        if shuffle_mode == "autologous":
            random_seq = "".join([str(i) for i in random.sample(seq, len(seq))])
        # Use transcriptome di-nucleotide frequencies for random sequences:
        elif shuffle_mode == "di-nuc":
            random_seq = generate_random_di_nuc_seq(len(seq))
        elif shuffle_mode == "transcriptome":
            random_seq = "".join([str(i) for i in random.sample(list_amount[num], len(seq))])
        randomized_sequences.append(random_seq)

    return randomized_sequences




###############################################################################################################################################
# CREATE SEQUENCE FILES FROM TRANSCRIPTOME
###############################################################################################################################################
def transcriptome_file(MANE_dict):
    subseqs = ["UTR3", "UTR5", "CDS", "transcript"]
    experiments = ["RNAcompete", "SELEX", "HT-SELEX"]
    transcriptome_3utr = "fimo_transcriptome_3utr.txt"
    transcriptome_5utr = "fimo_transcriptome_5utr.txt"
    transcriptome_cds = "fimo_transcriptome_cds.txt"
    transcriptome_full = "fimo_transcriptome_full.txt"

    transcriptome_files = [transcriptome_3utr,
                           transcriptome_5utr,
                           transcriptome_cds,
                           transcriptome_full]


    dir = os.path.join(os.getcwd(),"DATA", "FIMO_input", "sequences")
    if not os.path.isdir(dir):
        os.makedirs(dir)

    for i in range(len(transcriptome_files)):
        with open(os.path.join(dir, transcriptome_files[i]), "w") as f:
            for id in MANE_dict.keys():
                seq = MANE_dict[id][subseqs[i]]
                f.write(">" + id + " " + "transcriptome_" + subseqs[i] + "\n" + seq + "\n\n")



transcriptome_file(MANE_transcriptome)
#print(">>> DONE CREATING SEQUENCE FILES\n")


###############################################################################################################################################
# FUNCTIONS THAT WRITE GENE_ID:SEQUENCE DICTIONARIES INTO MOTIF-INPUT FILE
# "attract"-parameter is set to True when ATtRACT dict is entered (important for setting the file name)
###############################################################################################################################################


# Required form for input motifs:
# MEME version (MEME Suite 5.4.1)

# ALPHABET = RNA

# Background letter frequencies
# fasta-get-markov from MANE human sequences

# MOTIF Name/ID

def fimo_motif_input(dict_ppms, attract=True, rnacomp=True):
    if attract:
        if rnacomp:
            fimo_file = "fimo_attract_rnacomp.txt"
            exp = "RNAcompete"
        else:
            fimo_file = "fimo_attract_selex.txt"
            exp = "SELEX"
    else:
        fimo_file = "fimo_htselex.txt"
        exp = "HT-SELEX"

    dir = os.path.join(os.getcwd(), "DATA", "FIMO_input", "motifs")

    if not os.path.isdir(dir):
        os.makedirs(dir)

    open(os.path.join(dir, fimo_file), "a").close()

    with open(os.path.join(dir, fimo_file), "w") as f:
        f.write("MEME version 5.4.1\n\nALPHABET=ACGU\n")

        for i in range(len(list(dict_ppms))):

            id = list(dict_ppms)[i]
            matrix = dict_ppms[id]
            len_matrix = len(matrix)

            f.write("\nMOTIF "
                    + id
                    + " "
                    + exp
                    + "\nletter-probability matrix: alength= 4 "
                    + "w= "
                    + str(len(matrix) - 1)
                    + "\n")

            for k in range(len_matrix - 1):
                probabilities_list = matrix[k]
                f.write(" " + " ".join(map(str, probabilities_list)) + "\n")

                if k + 1 == len_matrix:
                    f.write("\n")
                    continue


fimo_motif_input(attract_ppms["RNAcomp"], attract=True, rnacomp=True)
fimo_motif_input(attract_ppms["SELEX"], attract=True, rnacomp=False)
fimo_motif_input(htselex_ppms, attract=False, rnacomp=False)

#print(">>> DONE CREATING MOTIF FILES\n")


'''###############################################################################################################################################
# SOME STATISTICS ON THE SEQUENCES; MEAN LENGTH; MIN AND MAX LENGTH
###############################################################################################################################################
import numpy as np
subseqs = ["UTR3", "UTR5", "CDS", "cDNA"]

seq_length = {}
lengths = []

for subseq in subseqs:
    seq_length[subseq] = {k:len(v) for (k,v) in MANE_transcriptome[subseq].items()}

stat_dict = {}

for subseq in subseqs:
    stat_dict[subseq] = {}
    stat_dict[subseq]["mean"] = np.mean(np.asarray(list(seq_length[subseq].values())))
    stat_dict[subseq]["min"] = np.min(np.asarray(list(seq_length[subseq].values())))
    stat_dict[subseq]["max"] = np.max(np.asarray(list(seq_length[subseq].values())))
    stat_dict[subseq]["median"] = np.median(np.asarray(list(seq_length[subseq].values())))
    stat_dict[subseq]["sum"] = np.sum(np.asarray(list(seq_length[subseq].values())))


exps = ["RNAcomp", "SELEX"]
ppm_mean_len = {}
ppm_lengths = []###############################################################################################################################################
# SOME STATISTICS ON THE SEQUENCES; MEAN LENGTH; MIN AND MAX LENGTH
###############################################################################################################################################
import numpy as np
subseqs = ["UTR3", "UTR5", "CDS", "cDNA"]

seq_length = {}
lengths = []

for subseq in subseqs:
    seq_length[subseq] = {k:len(v) for (k,v) in MANE_transcriptome[subseq].items()}

stat_dict = {}

for subseq in subseqs:
    stat_dict[subseq] = {}
    stat_dict[subseq]["mean"] = np.mean(np.asarray(list(seq_length[subseq].values())))
    stat_dict[subseq]["min"] = np.min(np.asarray(list(seq_length[subseq].values())))
    stat_dict[subseq]["max"] = np.max(np.asarray(list(seq_length[subseq].values())))
    stat_dict[subseq]["median"] = np.median(np.asarray(list(seq_length[subseq].values())))
    stat_dict[subseq]["sum"] = np.sum(np.asarray(list(seq_length[subseq].values())))


exps = ["RNAcomp", "SELEX"]
ppm_mean_len = {}
ppm_lengths = []
for exp in exps:
    for i in range(len(list(attract_ppms[exp].values()))):
        ppm_lengths.append(len(list(attract_ppms["RNAcomp"].values())[i]))
    ppm_mean_len[exp] = np.mean(ppm_lengths)

rnacomp_mean_len = ppm_mean_len["RNAcomp"]

tests = stat_dict["CDS"]["sum"]/rnacomp_mean_len
p_value = 1e-4/tests

for exp in exps:
    for i in range(len(list(attract_ppms[exp].values()))):
        ppm_lengths.append(len(list(attract_ppms["RNAcomp"].values())[i]))
    ppm_mean_len[exp] = np.mean(ppm_lengths)

rnacomp_mean_len = ppm_mean_len["RNAcomp"]

tests = stat_dict["CDS"]["sum"]/rnacomp_mean_len
p_value = 1e-4/tests
'''


