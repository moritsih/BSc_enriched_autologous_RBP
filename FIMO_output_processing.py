import os
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

###########################################################################
#CHANGE THESE SETTINGS TO RUN DIFFERENT ANALYSES
###########################################################################
pval_100 = False
pval_1000 = True
pval_10000 = False
transcriptome_background = True
individual_background = False
###########################################################################

if pval_100:
    diagram_title = "Enrichment of autologous binding via FIMO\nFull Transcriptome background\nP-value cutoff: 1e-2"
    data_path = os.path.join(os.getcwd(), "DATA", "FIMO_OUT", "background_transcriptome_pval1e-2")

if pval_1000:
    diagram_title = "Enrichment of autologous binding via FIMO\nFull Transcriptome background\nP-value cutoff: 1e-3"
    data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "background_transcriptome_pval1e-3")

if pval_10000:
    if transcriptome_background:
        diagram_title = "Enrichment of autologous binding via FIMO\nFull Transcriptome background\nP-value cutoff: 1e-3"
        data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "background_transcriptome_pval1e-4")

    if individual_background:
        diagram_title = "Enrichment of autologous binding via FIMO\nIndividual transcript parts as background\nP-value cutoff: 1e-4"
        data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "background_individual_pval1e-4")


directories = os.listdir(data_path)  # list of directiories containing output tsv files (& other files)

# used for loop
experiments = ["RNAcompete", "SELEX", "HT-SELEX"]
subsequences = ["UTR3", "UTR5", "CDS", "transcript"]

# prints amount of analyzed RBPs per experiment onto plot
rnacompete_ppms = list(attract_ppms["RNAcomp"].keys( ))
selex_ppms = list(attract_ppms["SELEX"].keys( ))
htselex_ppms = list(htselex_ppms.keys( ))
rbps = [rnacompete_ppms, selex_ppms, htselex_ppms]

# stores access to folders where the FIMO-output data lies; for retrieval inside loop
def store_filenames_for_retrieval():
    global data_path
    global directories

    matches_files = {}

    for dir in directories:  # loop through folders of different experiment/subsequence combinations
        exp = os.path.join(data_path, dir)  # define which experiment we're dealing with
        for file in os.listdir(exp):
            if ".tsv" in file:  # there are other types of files in the dir
                tsv_file_path = os.path.join(exp, file)
                matches_files[dir] = tsv_file_path

        print(f">>> FETCHED MATCHES FROM {dir}\n")

    return matches_files

matches_raw = store_filenames_for_retrieval( )


def sort_matches(matches_unsorted):
    sorter_exp = {"rnacomp": "RNAcompete", "selex": "SELEX", "htselex": "HT-SELEX"}
    sorter_subseq = {"UTR3": "UTR3", "UTR5": "UTR5", "CDS": "CDS", "full": "transcript"}

    matches_sorter = {}

    for exp_no, exp in sorter_exp.items( ):
        matches_sorter[exp] = {}

        for subseq_no, subseq in sorter_subseq.items( ):

            for dir in directories:
                if dir.startswith(exp_no) and dir.endswith(subseq_no):
                    matches_sorter[exp][subseq] = matches_unsorted[dir]

        print(f">>> {exp} IS IN RIGHT SHAPE NOW\n")

    return matches_sorter

matches_sorted = sort_matches(matches_raw)


def read_tsv_file(tsv_file_path):
    with open(tsv_file_path, "r") as f:
        _ = f.readline()
        content = f.read().split("\n")
        content = [x for x in content if x and not x.startswith("#")]  # removing bottom lines

        infos = []
        for line in content:
            line = line.split("\t")
            motif_id = line[0]
            seq_id = line[2]
            start = line[3]
            stop = line[4]
            infos.append([seq_id, motif_id, start, stop])

        return infos

def plot_num_matches(match_dict):
    global individual_background, transcriptome_background

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    if individual_background:
        fig.suptitle("Number of matches - background distributions by subsequences")
    if transcriptome_background:
        fig.suptitle("Number of matches - transcriptome background distribution")

    for i, exp in enumerate(matches_sorted.keys( )):
        subseq_matches = []
        seqs = list(matches_sorted[exp].keys( ))
        for seq in seqs:
            infos = read_tsv_file(matches_sorted[exp][seq])
            lengths = len(infos)
            subseq_matches.append(lengths)
        sns.set_style=("whitegrid")
        sns.barplot(ax=axes[i], x=seqs, y=subseq_matches)
        axes[i].set_title(exp)

    plt.show( )

#plot_num_matches(matches_sorted)



def pipeline_for_FIMO_analysis(matches_sorted_dict):
    global MANE_transcriptome
    global experiments
    global subsequences
    global attract_ppms, htselex_ppms
    global diagram_title

    SEED = 1234

    fig = plt.figure(figsize=[18.3 / 2.54, 11.0 / 2.54], constrained_layout=True, dpi=300)
    grid = fig.add_gridspec(1, 1)
    ax3 = plt.subplot(grid[0, 0])
    ax3.set_title(diagram_title)

    # great outer loop for going through files
    for i, exp in enumerate(experiments):

        for l, subseq in enumerate(subsequences):


            tsv_file_path = matches_sorted[exp][subseq]

            infos = read_tsv_file(tsv_file_path)

            add_sequence_length_to_infos(infos, subseq)

            infos, duplicated_matrices = group_by_motif_id_and_sequence_id(infos)

            #with tqdm(total=set_tqdm_counter_total()) as pbar:
            #    pbar.update(1)

            infos, \
            len_normalize, \
            consider_overlap, \
            normalize_by_num_of_matrices = calc_coverages_autol_len(infos,
                                                                    subseq,
                                                                    duplicated_matrices,
                                                                    SEED,
                                                                    MANE_transcriptome,
                                                                    len_normalize=False,
                                                                    consider_overlap=False,
                                                                    background_longer_than_autol_only=False,
                                                                    normalize_by_num_of_matrices=False)

            autologous_all_motifs = []
            background_all_motifs = []

            for motif in infos.keys():
                mean_cov_per_motif, std_cov_per_motif = get_mean_std(infos[motif])
                autologous, background = calc_z_scores(infos[motif],
                                                       mean_cov_per_motif,
                                                       std_cov_per_motif)

                autologous_all_motifs.append(autologous)
                background_all_motifs.append(background)

            p_val = binary_vs_averaged(autologous_all_motifs, background_all_motifs, ranked="no", side="higher")

            number_of_motifs_used = count_amount_of_motifs_per_experiment(attract_ppms, htselex_ppms)

            set_plotting_details()

            plot_analysis_results(ax3,
                                  autologous_all_motifs,
                                  background_all_motifs,
                                  p_val,
                                  number_of_motifs_used,
                                  i, # see for-loop above: i helps put labels on the plot
                                  l,
                                  subseq)



        print(">>> DONE WITH COVERAGE CALCULATIONS\n")
        print(">>> Z-SCORES HAVE BEEN CALCULATED\n")
        print(">>> PLOTTING THE DATA\n")

    plt.grid()
    plt.show()


def set_tqdm_counter_total(match_list):
    return len(match_list)


def set_plotting_details():

    size_tiny = 3
    size_small = 5
    size_medium = 9
    size_big = 10
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["font.family"] = "Calibri"
    plt.rcParams["lines.markersize"] = 4
    plt.rc("font", size=size_tiny)          # controls default text sizes
    plt.rc("axes", titlesize=size_tiny)     # fontsize of the axes title
    plt.rc("axes", labelsize=size_tiny)     # fontsize of the x and y labels
    plt.rc("xtick", labelsize=size_small)    # fontsize of the tick labels
    plt.rc("ytick", labelsize=size_small)    # fontsize of the tick labels
    plt.rc("legend", fontsize=size_tiny)     # legend fontsize
    plt.rc("figure", titlesize=size_big)     # fontsize of the figure title
    mpl.rcParams["legend.markerscale"] = 1.0


def count_amount_of_motifs_per_experiment(attract, htselex):

    number_of_motifs_used = []
    for k in attract.keys():
        number_of_motifs_used.append(len(attract[k]))

    number_of_motifs_used.append(len(htselex))

    return number_of_motifs_used



def add_sequence_length_to_infos(infos, subseq):
    global MANE_transcriptome

    for i in range(len(infos)):
        seq_id = infos[i][0]
        seq_len = len(MANE_transcriptome[seq_id][subseq])
        infos[i].append(seq_len)



def group_by_motif_id_and_sequence_id(matches_by_motif, merge_duplicate_motifs=False):
    motif_subseq = {}
    duplicated_matrices = {}

    for match in matches_by_motif:
        duplication_counter = 1 # counts amount of different matrices per motif (if there are duplicates)
        seq_id = match[0]
        motif_id = match[1]
        if merge_duplicate_motifs: # Some RNAcompete RBPs had multiple matrices; Merge matches to 1 or keep separate?
            if "_" in motif_id:
                num_of_duplicates = int(motif_id[-1]) + 1 # if matrix has suffix _4, then it'll store "5"
                motif_id = motif_id[:-2] #removes _1 suffix; enables "merging" of matches
                duplicated_matrices[motif_id] = num_of_duplicates

        if motif_id not in motif_subseq:
            motif_subseq[motif_id] = {}

        if seq_id not in motif_subseq[motif_id]:
            motif_subseq[motif_id][seq_id] = []

        motif_subseq[motif_id][seq_id].append(match)
    return motif_subseq, duplicated_matrices




def calc_coverages_autol_len(matches_by_subseq,
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
# ppms gives the number of matching motifs to the plot (N = )
    # coverages: holds motifs and the sequences each motif matched with. For every sequence, there's a coverage



def get_mean_std(matches_by_sequences):

    mean_cov = []
    std_cov = []

    for box in matches_by_sequences.values():
        #box = all sequences a motif matched with
        mean_cov.append(box[2])
        std_cov.append(box[2])

    mean_cov = np.nanmean(mean_cov)
    std_cov = np.nanstd(std_cov)

    return mean_cov, std_cov



def calc_z_scores(coverages_by_motif, mean_cov, std_cov):

    mean = mean_cov
    std = std_cov
    background = []
    autologous_match_occurred = False


    for info in coverages_by_motif.values():

        seq_id = info[0]
        motif_id = info[1]
        z_val = (info[2] - mean) / std
        background.append(z_val)

        if seq_id == motif_id:
            autologous = z_val
            autologous_match_occurred = True

    if not autologous_match_occurred:
        autologous = ((0 - mean) / std)

    return autologous, background


def pie_plot():
    subsequences = ["autologous UTR3", "autologous UTR5", "autologous CDS", "autologous transcript"]

    fig = plt.figure()

    colors = ['b', 'g', 'r', 'c']
    global experiments, final_distributions
    for i, exp in enumerate(experiments):
        x = 131
        ax = fig.add_subplot(x+i)
        data = []
        for subseq in subsequences:
            data.append(len(final_distributions[exp][subseq]))
        ax.set_title(exp)
        ax.pie(data, labels=subsequences,
               colors=colors,
               textprops={'fontsize': 8},
               autopct='%.1f%%')
    #plt.tight_layout()
    plt.show()

#pie_plot()

def z_score_plot():
    global experiments, final_distributions
    subsequences = ["UTR3", "UTR5", "CDS", "transcript",
                    "autologous UTR3", "autologous UTR5", "autologous CDS", "autologous transcript"]
    fig = plt.subplot(len(experiments), )

    grid = fig.add_gridspec(len(experiments), len(subsequences))


    for i,exp in enumerate(experiments):
        for l,subseq in enumerate(subsequences):

            ax = fig.add_subplot(grid[i, l])
            scores = final_distributions[exp][subseq]
            x = np.linspace(-6, 6, num=len(final_distributions[exp][subseq]))
            ax.scatter(x = x, y = scores, label=[exp, subseq])

    plt.legend()
    plt.show()

#z_score_plot()

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
                    label="autol.(z-score in resp. transcriptome)")


        ax3.plot([i-0.25,i-0.05],
                 [np.mean(autologous)]*2,
                 color=point_colors[l],
                 zorder=2,
                 label="(autologous)(mean)",
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
                 label="(merged histogram)",
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
             f"$p_{subseq}$ = {pvalue}",
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