import os
import sys
from matplotlib import pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy import stats
import random
from FIMO_input_processing import MANE_transcriptome, attract_ppms, htselex_ppms
MANE_transcriptome = MANE_transcriptome
from additional_code.load_histograms import vertical_hist
import seaborn as sns
from tqdm import tqdm
import mmap

###########################################################################
#CHANGE THESE SETTINGS TO RUN DIFFERENT ANALYSES
###########################################################################
pval_cutoff = sys.argv[1]
###########################################################################

if pval_cutoff == "5e-2":
    diagram_title = "Enrichment of autologous binding via FIMO\nFull Transcriptome background\nP-value cutoff: 0.05"
    data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "pval5e-2")
    figure_name = "no_multiple_motifs_pval5e-2"

if pval_cutoff == "1e-2":
    diagram_title = "Enrichment of autologous binding via FIMO\nFull Transcriptome background\nP-value cutoff: 1e-2"
    data_path = os.path.join(os.getcwd(), "DATA", "FIMO_OUT", "pval1e-2")
    figure_name = "no_multiple_motifs_pval1e-2"

if pval_cutoff == "1e-3":
    diagram_title = "Enrichment of autologous binding via FIMO\nFull Transcriptome background\nP-value cutoff: 1e-3"
    data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "pval1e-3")
    figure_name = "no_multiple_motifs_pval1e-3"

if pval_cutoff == "1e-4":
    diagram_title = "Enrichment of autologous binding via FIMO\nFull Transcriptome background\nP-value cutoff: 1e-4"
    data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "pval1e-4")
    figure_name = "no_multiple_motifs_pval1e-4"

files = os.listdir(data_path)  # list of directiories containing output tsv files (& other files)
plot_target_path = os.path.join(os.getcwd( ), "DATA", "result_plots")

# used for loop
experiments = ["RNAcompete", "SELEX", "HT-SELEX"]
subsequences = ["UTR3", "UTR5", "CDS", "transcript"]

# prints amount of analyzed RBPs per experiment onto plot
rnacompete_ppms = list(attract_ppms["RNAcomp"].keys( ))
selex_ppms = list(attract_ppms["SELEX"].keys( ))
htselex_ppms = list(htselex_ppms.keys( ))
rbps = [rnacompete_ppms, selex_ppms, htselex_ppms]





# stores access to folders where the FIMO-output data lies; for retrieval inside loop
def store_filenames_for_retrieval(data_path, files):

    matches_files = {}

    for file in files:  # loop through files with different experiment/subsequence combinations
        tsv_file_path = os.path.join(data_path, file)
        matches_files[file] = tsv_file_path

        print(f">>> FETCHED MATCHES FROM {file}\n")

    return matches_files

matches_raw = store_filenames_for_retrieval(data_path, files)


def sort_matches(matches_unsorted):
    sorter_exp = {"rnacomp": "RNAcompete", "selex": "SELEX", "htselex": "HT-SELEX"}
    sorter_subseq = {"UTR3": "UTR3", "UTR5": "UTR5", "CDS": "CDS", "full": "transcript"}

    matches_sorter = {}

    for exp_no, exp in sorter_exp.items( ):
        matches_sorter[exp] = {}

        for subseq_no, subseq in sorter_subseq.items( ):

            for file in files:
                name = str(file)[:-4]
                if name.startswith(exp_no) and name.endswith(subseq_no):
                    matches_sorter[exp][subseq] = matches_unsorted[file]

    return matches_sorter

matches_sorted = sort_matches(matches_raw)








def pipeline_for_FIMO_analysis(matches_sorted_dict):
    global MANE_transcriptome
    global experiments
    global subsequences
    global attract_ppms, htselex_ppms
    global diagram_title
    global plot_target_path

    SEED = 1234

    fig = plt.figure(figsize=[18.3 / 2.54, 11.0 / 2.54], constrained_layout=True, dpi=300)
    grid = fig.add_gridspec(1, 1)
    ax3 = plt.subplot(grid[0, 0])
    ax3.set_title(diagram_title)

    set_plotting_details( )

    # great outer loop for going through files
    for i, exp in enumerate(experiments):
        print(f"STARTING ON {exp} ...")

        for l, subseq in enumerate(subsequences):
            print(f">>> STEPPING INTO {subseq} ...")

            tsv_file_path = matches_sorted[exp][subseq]

            print(">>> EXTRACTING MATCHES FROM FILES ...")

            num_of_lines = get_number_of_lines(tsv_file_path)

            print(">>> GOT NUMBER OF MATCHES ...")

            info_generator = file_info_generator(tsv_file_path)

            array_per_motif_seq_combination = {}

            print(">>> LARGE FILES ARE BEING PROCESSED ...")
            for infos in tqdm(info_generator, total=num_of_lines):
            #for infos in info_generator:

                seq_bucket, start, stop = create_array_of_seq_length(infos,
                                                                     array_per_motif_seq_combination,
                                                                     subseq,
                                                                     MANE_transcriptome)

                flip_zeros_to_ones(seq_bucket, start, stop, consider_overlap=False)

            print(">>> COVERAGES ARE BEING COMPUTED ...")
            infos = calculate_coverage_per_seq(array_per_motif_seq_combination,
                                               MANE_transcriptome,
                                               subseq,
                                               length_normalization=False)

            autologous_all_motifs = []
            background_all_motifs = []

            print(">>> MEAN+STDS ARE BEING CALCULATED AND Z-SCORES ARE BEING COMPUTED ...")
            for motif in infos.keys():
                mean_cov_per_motif, std_cov_per_motif = get_mean_std(infos[motif])
                autologous, background = calc_z_scores(infos[motif],
                                                       motif,
                                                       mean_cov_per_motif,
                                                       std_cov_per_motif)

                autologous_all_motifs.append(autologous)
                background_all_motifs.append(background)

            print(">>> P-VALUES ARE BEING CALCULATED ...")
            p_val = binary_vs_averaged(autologous_all_motifs, background_all_motifs, ranked="no", side="higher")

            number_of_motifs_used = count_amount_of_motifs_per_experiment(attract_ppms, htselex_ppms)

            print(">>> PUTTING INTO PLOT ...")
            plot_analysis_results(ax3,
                                  autologous_all_motifs,
                                  background_all_motifs,
                                  p_val,
                                  number_of_motifs_used,
                                  i, # see for-loop above: i helps put labels on the plot
                                  l,
                                  subseq)


    plt.grid()
    plt.savefig(os.path.join(plot_target_path,figure_name))
    plt.show()


def get_number_of_lines(tsv_file_path):
    with open(tsv_file_path, "r") as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            total_lines = 0
            while mm.readline():
                total_lines += 1
    return total_lines


def file_info_generator(tsv_file_path):

    with open(tsv_file_path, "r") as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            _ = mm.readline()
            while True:
                line = mm.readline()
                if not line:
                    break
                line = line.decode('UTF-8')
                if not line.startswith('#'):
                    line = line.split("\t")
                    motif_id = line[0]
                    seq_id = line[2]
                    start = line[3]
                    stop = line[4]
                    infos = [seq_id, motif_id, start, stop]

                    yield infos



def create_array_of_seq_length(infos, array_sorter, subseq, MANE_transcriptome):

    seq_id = infos[0]
    matrix_id = infos[1]
    start = infos[2]
    stop = infos[3]
    seq_len = len(MANE_transcriptome[seq_id][subseq])
    if matrix_id not in array_sorter:
        array_sorter[matrix_id] = {}

    matrix_bucket = array_sorter[matrix_id]
    if seq_id not in matrix_bucket:
        matrix_bucket[seq_id] = np.zeros(seq_len)

    return array_sorter[matrix_id][seq_id], int(start), int(stop)


def flip_zeros_to_ones(zero_array, start, stop, consider_overlap=False):

    if not consider_overlap:
        zero_array[start:stop] = 1
    else:
        zero_array[start:stop] =+ 1


def calculate_coverage_per_seq(cov_arrays, MANE_transcriptome, subseq, length_normalization=False):
    for matrix in cov_arrays.keys( ):
        matrix_bucket = cov_arrays[matrix]

        for seq in matrix_bucket.keys( ):
            seq_len = len(MANE_transcriptome[seq][subseq])

            if length_normalization:
                matrix_bucket[seq] = np.mean(matrix_bucket[seq]) / seq_len
            else:
                matrix_bucket[seq] = np.mean(matrix_bucket[seq])

    return cov_arrays



def set_tqdm_counter_total(match_list):
    return len(match_list)


def set_plotting_details():

    size_tiny = 3
    size_small = 5
    size_medium = 9
    size_big = 12
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["font.family"] = "Calibri"
    plt.rcParams["lines.markersize"] = 4
    plt.rc("font", size=size_small)          # controls default text sizes
    plt.rc("axes", titlesize=size_small)     # fontsize of the axes title
    plt.rc("axes", labelsize=size_small)     # fontsize of the x and y labels
    plt.rc("xtick", labelsize=size_small)    # fontsize of the tick labels
    plt.rc("ytick", labelsize=size_small)    # fontsize of the tick labels
    plt.rc("legend", fontsize=size_small)     # legend fontsize
    plt.rc("figure", titlesize=size_big)     # fontsize of the figure title
    mpl.rcParams["legend.markerscale"] = 1.0




def count_amount_of_motifs_per_experiment(attract, htselex):

    number_of_motifs_used = []
    for k in attract.keys():
        number_of_motifs_used.append(len(attract[k]))

    number_of_motifs_used.append(len(htselex))

    return number_of_motifs_used



def add_sequence_length_to_infos(infos, subseq, MANE_transcriptome):

    seq_id = infos[0]
    seq_len = len(MANE_transcriptome[seq_id][subseq])
    infos.append(seq_len)




def calc_coverages(matches_by_subseq,
                   subseq,
                   duplicated_matrices,
                   SEED,
                   MANE_transcriptome,
                   len_normalize=False,
                   consider_overlap=False,
                   background_longer_than_autol_only=False,
                   normalize_by_num_of_matrices=False):

    np.random.seed(SEED)
    coverages = {}

    for motif in matches_by_subseq:
        coverages[motif] = {}

        for seq in matches_by_subseq[motif].values(): #all the matches a motif has with a sequence

            autologous_seq_len = len(MANE_transcriptome[motif][subseq]) # get the length of the autologous sequence
            seq_len = seq[0][4]
            seq_id = seq[0][0]
            motif_id = seq[0][1]

            if background_longer_than_autol_only: # only include background if the matched sequence is of same length
                                                  # or longer than autologous sequence

                if seq_len >= autologous_seq_len:

                    for idx in range(len(seq)): # how many matches are there with a single sequence? idx will tell you

                        data = seq[idx]
                        start = int(data[2])
                        stop = int(data[3])


                        if idx == 0:
                            array = np.zeros(seq_len)
                            array[start - 1:stop] = 1 # numpy 1:5 means= 2,3,4,5;

                        else:

                            if consider_overlap:
                                array[start - 1:stop] += 1 # overlaps get added and count more

                            else:
                                array[start - 1:stop] = 1


                    start_idx_range = seq_len-autologous_seq_len

                    if start_idx_range == 0:
                        cov = [seq_id, motif_id, np.mean(array)]
                        continue


                    start_idx = np.random.randint(0, start_idx_range)
                    stop_idx = start_idx + autologous_seq_len

                    array = array[start_idx:stop_idx]

                    cov = [seq_id, motif_id, np.mean(array)]

                else:
                    continue

                if normalize_by_num_of_matrices and len_normalize:
                    if motif in duplicated_matrices:
                        num_of_duplicates = duplicated_matrices[motif]
                    else:
                        num_of_duplicates = 1
                    cov[2] = cov[2] / (num_of_duplicates * seq_len)
                    coverages[motif][seq_id] = cov
                    continue

                if normalize_by_num_of_matrices:
                    if motif in duplicated_matrices:
                        num_of_duplicates = duplicated_matrices[motif]
                    else:
                        num_of_duplicates = 1
                    cov[2] = cov[2] / num_of_duplicates
                    coverages[motif][seq_id] = cov
                    continue

                if len_normalize:
                    cov[2] = cov[2] / seq_len
                    coverages[motif][seq_id] = cov
                    continue

                else:
                    coverages[motif][seq_id] = cov



            else:
                for idx in range(len(seq)):

                    data = seq[idx]
                    seq_id = data[0]
                    motif_id = data[1]
                    start = int(data[2])
                    stop = int(data[3])
                    seq_len = data[4]

                    if idx == 0:
                        array = np.zeros(seq_len)
                        array[start - 1:stop] = 1  # numpy 1:5 means= 2,3,4,5;

                    else:

                        if consider_overlap:
                            array[start - 1:stop] += 1  # overlaps get added and count more

                        else:
                            array[start - 1:stop] = 1

                cov = [seq_id, motif_id, np.mean(array), seq_len]

                if normalize_by_num_of_matrices and len_normalize:
                    if motif in duplicated_matrices:
                        num_of_duplicates = duplicated_matrices[motif]
                    else:
                        num_of_duplicates = 1
                    cov[2] = cov[2] / (num_of_duplicates * seq_len)
                    coverages[motif][seq_id] = cov
                    continue

                if normalize_by_num_of_matrices:
                    if motif in duplicated_matrices:
                        num_of_duplicates = duplicated_matrices[motif]
                    else:
                        num_of_duplicates = 1
                    cov[2] = cov[2] / num_of_duplicates
                    coverages[motif][seq_id] = cov
                    continue

                if len_normalize:
                    cov[2] = cov[2] / seq_len
                    coverages[motif][seq_id] = cov
                    continue

                else:
                    coverages[motif][seq_id] = cov


    return coverages, len_normalize, consider_overlap, normalize_by_num_of_matrices
# coverages = data; len_norm and consider_overlap are for specifying the mode of analysis during plotting;
# coverages: holds motifs and the sequences each motif matched with. For every sequence, there's a coverage



def get_mean_std(matches_by_sequences):

    mean_cov = []
    std_cov = []

    for box in matches_by_sequences.values():
        mean_cov.append(box)
        std_cov.append(box)

    mean_cov = np.mean(mean_cov)
    std_cov = np.std(std_cov)

    return mean_cov, std_cov



def calc_z_scores(coverages_by_motif, motif_id, mean_cov, std_cov):

    mean = mean_cov
    std = std_cov
    background = []
    autologous_match_occurred = False


    for seq_id, cov in coverages_by_motif.items():

        z_val = (cov - mean) / std
        background.append(z_val)

        if seq_id == motif_id:
            autologous = z_val
            autologous_match_occurred = True

    if not autologous_match_occurred:
        autologous = ((0 - mean) / std)

    return autologous, background


my_EXPERIMENTS = ["RNAcompete", "SELEX", "HT-SELEX"]


# Function for getting p-value of autologous transcript:
def binary_vs_averaged(bi_ls, bi_ls_ls, ranked="no", side="higher"):
    if ranked == "yes":
        singles = []
        dists = []
        for i,element in enumerate(bi_ls_ls):
            ranks = stats.rankdata(element + bi_ls[i], method="average")/(len(element)+1)
            singles.append(ranks[-1])
            dists.append(ranks[:-1])
        bi_ls = singles
        bi_ls_ls = dists

    bi_total = sum(bi_ls)
    dist = []

    N=10**5
    hist = []
    for _ in range(N):
        total = 0
        for item in bi_ls_ls:
            chosen_one = random.choice(item)
            total += chosen_one
            hist.append(chosen_one)
        dist.append(total)

    lower, equal, higher = 0,0,0
    for element in dist:
        if element > bi_total:
            higher += 1
        if element < bi_total:
            lower += 1
        if element == bi_total:
            equal += 1

    if higher > lower:
        p_1 = lower/N + equal/N/2
        p_2 = p_1 *2

    if higher < lower:
        p_1 = higher/N + equal/N/2
        p_2 = p_1 *2

    if side == "lower":
        p_1 = lower/N + equal/N/2
        return p_1#, hist

    if side == "higher":
        p_1 = higher/N + equal/N/2
        return p_1#, hist

    return p_2#, hist


# With this function, a combined plot for the different experiments is created and the p values
# are calculated. The plot is then shown on screen.
def plot_analysis_results(ax3, autologous, background, pvalue, ppms, i, l, subseq):

    # As different subsequences need slighly different settings, I use "l" to go through these lists
    autologous_points_horizontal_step = [0.05, 0.1, 0.15, 0.2]
    point_colors = ["red", "orange", "lightblue", "blue"]
    p_value_position = [-7, -7.5, -8, -8.5]


    ##########################################################
    # PLOT LEFT SIDE of individual diagrams (autologous part):
    ##########################################################
    if i == 0:



        ax3.scatter([i-autologous_points_horizontal_step[l]] * len(autologous),
                    autologous,
                    linewidths=0,
                    color=point_colors[l],
                    alpha=0.2,
                    zorder=3,
                    label=f"autol. {subseq} (z-score in resp. transcriptome)")


        ax3.plot([i-0.25,i-0.05],
                 [np.mean(autologous)]*2,
                 color=point_colors[l],
                 zorder=2,
                 label=f"(autologous) {subseq} (mean)",
                 alpha=0.4)


    else:
        ax3.scatter([i - autologous_points_horizontal_step[l]] * len(autologous),
                    autologous,
                    linewidths=0,
                    color=point_colors[l],
                    alpha=0.2,
                    zorder=3)


        ax3.plot([i - 0.25, i - 0.05],
                 [np.mean(autologous)] * 2,
                 color=point_colors[l],
                 zorder=2,
                 alpha=0.4)


    hist = vertical_hist(background,
                         i+0.0,
                         0.6,
                         borders=[-6,6],
                         bin_number="standard",
                         window_average=5)


    vxs, vys, mean = hist["x_smooth"],hist["y_smooth"],hist["mean"]


    ##############################################################
    # PLOT RIGHT SIDE of individual diagrams (transcriptome part):
    ##############################################################
    if i == 0:

        ax3.plot(vxs,
                 vys,
                 color=point_colors[l],
                 zorder=2,
                 label=f"{subseq} (merged histogram)",
                 alpha=0.5,
                 linewidth=1)

        ax3.plot([i+0.05,i+0.25],
                 [mean]*2,
                 color=point_colors[l],
                 zorder=2,
                 alpha=0.5)

    else:

        ax3.plot(vxs,
                 vys,
                 color=point_colors[l],
                 zorder=2,
                 alpha=0.5,
                 linewidth=1)

        ax3.plot([i+0.05,i+0.25],
                 [mean]*2,
                 color=point_colors[l],
                 zorder=2,
                 alpha=0.5)


    #Write number of RBPs and p values for different experiments:
    ax3.text(i,
             -6.5,
             "N = %.0f"%ppms[i],
             ha="center")


    ax3.text(i,
             p_value_position[l],
             f"$p-{subseq}$ = {pvalue}",
             ha="center",
             color=point_colors[l])


    ax3.set_ylim([-6,6])

    ticklist = [k for k in range(len(my_EXPERIMENTS))]

    ax3.xaxis.tick_top()

    ax3.set_xticks(ticklist)

    ax3.set_xticklabels(my_EXPERIMENTS)

    ax3.legend(ncol=2,
               edgecolor="grey",
               framealpha = 0.8,
               facecolor = "white",
               frameon=True,
               loc=4)

    ax3.set_ylabel("motif coverage (nt/nt), z-score")



pipeline_for_FIMO_analysis(matches_sorted)