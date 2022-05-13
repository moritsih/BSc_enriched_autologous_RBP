from additional_code.read_MANE import *
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter

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
    gene_ids_with_duplicates = []
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

                elif gene_id in gene_ids_with_duplicates:
                    counter_for_matrix_duplicates = gene_ids_with_duplicates.count(gene_id)
                    gene_ids_with_duplicates.append(gene_id)
                    gene_id = gene_id + "_" + str(counter_for_matrix_duplicates)
                    dict_htselex_ppms[gene_id] = dict_curr_PPM
                    continue

                else:
                    gene_ids_with_duplicates.append(gene_id)
                    dict_htselex_ppms[gene_id] = dict_curr_PPM



    return dict_htselex_ppms

data_path = os.path.abspath("DATA")

attract_ppms = extract_attract_ppm(os.path.join(data_path, "ATtRACT_ppm.txt"))
htselex_ppms = extract_selex_ppm(os.path.join(data_path, "SELEX2020_ppm.txt"))

print(">>> DONE EXTRACTING PPMS FROM ATTRACT AND HT-SELEX\n")

###############################################################################################################################################
# FUNCTION GETTING NESTED DICT FROM ATtRACT DATABSE:
# EXPERIMENT:GENE_ID and GENE_ID:SEQUENCE
# + ONLY INCLUDES GENES THAT ARE AVAILABLE IN MANE DATABASE

# IMPORTANT VARIABLE: attract_ppms
###############################################################################################################################################


def read_db_attract(f_organism, f_experiment, f_score=True, f_mutated=False):

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
            #Gene_name = line[0]
            Gene_id = line[1]
            Mutated = line[2]
            Organism = line[3]
            #Motif = line[4]
            #Len = line[5]
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
            if f_score: # this filter is important for the inclusion of multiple matrices per protein;
                if "**" not in Score:
                    continue

            # Optional filters:

            if f_mutated:
                if "yes" in Mutated:
                    continue

            if Gene_id not in dict_MANE:
                continue

            if Gene_id in gene_ids_with_duplicates:
                if f_experiment == "RNAcompete": # RNAcompete entries only have 1 matrix/protein that has score with "**"
                    # therefore all but 1 entry are filtered out, allowing for this way of including multiple matrices/protein
                    counter_for_matrix_duplicates = gene_ids_with_duplicates.count(Gene_id)
                    gene_ids_with_duplicates.append(Gene_id)
                    Gene_id = Gene_id+"_"+str(counter_for_matrix_duplicates)
                    dict_attract_db[Gene_id] = Matrix_id
                    continue

                if f_experiment == "SELEX":
                    if Matrix_id in dict_attract_db.values():
                        continue
                    else:
                        counter_for_matrix_duplicates = gene_ids_with_duplicates.count(Gene_id)
                        gene_ids_with_duplicates.append(Gene_id)
                        Gene_id = Gene_id + "_" + str(counter_for_matrix_duplicates)
                        dict_attract_db[Gene_id] = Matrix_id
                        continue

            gene_ids_with_duplicates.append(Gene_id)

            dict_attract_db[Gene_id] = Matrix_id

        attract_ppms = {}
        for gene_id, matrix_id in dict_attract_db.items():
            attract_ppms[gene_id] = dict_ppms[matrix_id]

    return attract_ppms


attract_ppms = {}

attract_ppms["RNAcomp"] = read_db_attract("Homo_sapiens", "RNAcompete", f_score=True, f_mutated=False)
attract_ppms["SELEX"] = read_db_attract("Homo_sapiens", "SELEX", f_score=True, f_mutated=False)


RNAcomp_ppm_num = len(attract_ppms["RNAcomp"].keys())
SELEX_ppm_num = len(attract_ppms["SELEX"].keys())
htselex_ppm_num = len(htselex_ppms.keys())



print(">>> DONE CREATING PPM DICTIONARIES FOR FURTHER PROCESSING\n")
print(f">>> NUM OF ATTRACT RNAcomp MOTIFS (pre-filter): {RNAcomp_ppm_num}\n")
print(f">>> NUM OF ATTRACT SELEX MOTIFS (pre-filter): {SELEX_ppm_num}\n")
print(f">>> NUM OF HT-SELEX MOTIFS (pre-filter): {htselex_ppm_num}\n")


###############################################################################################################################################
# PLOTTING LENGTH DISTRIBUTION OF MOTIFS PER EXPERIMENT
###############################################################################################################################################

ppms = {}
ppms["RNAcompete"] = attract_ppms["RNAcomp"]
ppms["SELEX"] = attract_ppms["SELEX"]
ppms["HT-SELEX"] = htselex_ppms

def get_list_of_motif_ids(experiment):

    matrix_ids = [x if "_" not in x else x[:-2] for x in list(experiment.keys())]
    return matrix_ids

def count_number_of_matrices_per_protein(matrix_id_list):
    return Counter(matrix_id_list)

def plot_matrices_per_protein(ppms):

    fig, ax = plt.subplots(1, 3, figsize=(10,6))
    fig.suptitle("Number of matrices per protein", size=23)
    colors = ["red", "orange", "lightblue"]

    for i, exp in enumerate(ppms.keys()):
        counted_matrices = count_number_of_matrices_per_protein(
            get_list_of_motif_ids(
                ppms[exp]
            )
        )
        bins, counts = extract_plotable_data(counted_matrices.values())

        ax[i].bar(x=bins, height=counts, label=exp, color=colors[i])
        ax[i].set_xticks(np.arange(0, max(bins) + 1, 1))
        ax[1].set_xlabel("Amount of matrices", labelpad=10.0, size=17)
        ax[0].set_ylabel("Number of proteins", labelpad=15.0, size=17)
        ax[i].set_yticks(np.arange(0, max(counts), round(max(counts)/2)))
        ax[i].legend(loc="upper right", prop={'size': 10})

    plt.tight_layout()

    fig.savefig("DATA/result_plots/number_of_matrices_per_protein.png")
    #plt.show()


def extract_plotable_data(counted_matrices):
    counter = Counter(counted_matrices)
    bins = np.arange(min(counter.keys()), max(counter.keys())+1, 1)
    counts = [0 if x not in counter.keys() else counter[x] for x in bins]
    return bins, counts


plot_matrices_per_protein(ppms)

def extract_list_of_motif_lengths(ppms):
    dict_for_lists = {}
    for exp, motifs in ppms.items( ):
        dict_for_lists[exp] = sorted(map(len, motifs.values()))
    return dict_for_lists


def count_occurrences_of_motif_length(dict_motif_lengths):
    dict_for_counters = {}
    for exp, motif_lengths in dict_motif_lengths.items( ):
        dict_for_counters[exp] = Counter(motif_lengths)
    return dict_for_counters


def make_counter_into_barplot_data(counter_dict):
    ele_list = [int(x) for x in list(counter_dict.keys( ))]
    min_ele = np.min(ele_list)
    max_ele = np.max(ele_list)
    keys = np.arange(min_ele, max_ele + 1, 1)
    values = [0 if x not in counter_dict.keys( ) else counter_dict[x] for x in keys]

    return keys, values


def plot_motif_lengths(counter_dict):
    fig, ax = plt.subplots(3, 1, figsize=(10, 6))
    fig.suptitle("Motif lengths by experiment", size=25)
    colors = ["red", "orange", "lightblue"]


    for i, exp in enumerate(counter_dict.keys( )):

        x, y = make_counter_into_barplot_data(counter_dict[exp])

        ax[i].bar(x, y, width = 0.5, label=exp, color=colors[i], alpha=0.5, linewidth=0.5)
        ax[2].set_xlabel("Distribution of matrix lengths", labelpad=10.0, size=17)
        ax[1].set_ylabel("Matrix occurrences of given length", labelpad=15.0, size=17)
        ax[i].set_xticks(np.arange(0, max(x)+3, 1))
        ax[i].set_yticks(np.arange(0, max(y), round(max(y)/4)))
        ax[i].legend(loc="upper right", prop={'size': 10})

    plt.tight_layout()
    fig.savefig("DATA/result_plots/motif_length_distributions.png")
    #plt.show( )

plot_motif_lengths(count_occurrences_of_motif_length(extract_list_of_motif_lengths(ppms)))


def plot_num_of_motifs_per_protein(ppms):
    from collections import Counter
    fig, ax = plt.subplots(1,3, figsize=(18,5))
    for i, ppm in enumerate(ppms):

        motif_list = []
        for motif in ppm.keys():
            if "_" in motif:
                motif = str(motif)[:-2]
            motif_list.append(motif)

        matrices_per_motif = Counter(motif_list)

        matrices = matrices_per_motif.keys()
        count = matrices_per_motif.values()


        ax[i].hist(count, matrices)

    plt.show()

#plot_num_of_motifs_per_protein([attract_ppms["RNAcomp"], attract_ppms["SELEX"], htselex_ppms])





###############################################################################################################################################
# FILTER TRANSCRIPTOME DATABASE FOR UTR3, UTR5, CDS AND FULL TRANSCRIPT
###############################################################################################################################################


subseqs = ["UTR3", "UTR5", "CDS", "cDNA"]
MANE_transcriptome = {}

for id, content in dict_MANE.items():
    MANE_transcriptome[id] = {}
    for subseq in subseqs:
        MANE_transcriptome[id][subseq] = dict_MANE[id][subseq]

print(">>> DONE SPLITTING MANE SEQUENCES INTO SUBSEQUENCES\n")


def rename_MANE_seqs_duplicated_motifs(MANE_transcriptome, ppms):
    for ppm in ppms:
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
MANE_transcriptome = rename_MANE_seqs_duplicated_motifs(MANE_transcriptome, attract_ppms["SELEX"])
MANE_transcriptome = rename_MANE_seqs_duplicated_motifs(MANE_transcriptome, htselex_ppms)


###############################################################################################################################################
# FILTER UTRs AND CDS FOR EMPTY/SHORT SEQUENCES;
# DETERMINE A CUTOFF FOR SEQUENCE LENGTH 20
###############################################################################################################################################
def len_cutoff_subseq(cutoff=20):

    for id in MANE_transcriptome.keys():
        for subseq in MANE_transcriptome[id].keys():
            if len(MANE_transcriptome[id][subseq]) <= cutoff:
                MANE_transcriptome[id][subseq] = ""

    print(">>> DONE FILTERING FOR SUBSEQUENCE LENGTH\n")

len_cutoff_subseq(cutoff=0) #cutoff is not necessary after all since FIMO does this for us



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
            for id in [x for x in MANE_dict.keys() if "_" not in x]:
                seq = MANE_dict[id][subseqs[i]]
                f.write(">" + id + " " + subseqs[i] + "\n" + seq + "\n\n")



transcriptome_file(MANE_transcriptome)
print(">>> DONE CREATING SEQUENCE FILES\n")


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

print(">>> DONE CREATING MOTIF FILES\n")
