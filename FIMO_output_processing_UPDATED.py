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


###########################################################################
#CHANGE THESE SETTINGS TO RUN DIFFERENT ANALYSES
###########################################################################
pval_100 = False
pval_1000 = False
pval_10000 = True
transcriptome_background = True
individual_background = False
###########################################################################
if pval_100:
    diagram_title = "Enrichment of autologous binding via FIMO\nFull Transcriptome background of first order markov chain\nP-value cutoff: 1e-2"
    data_path = os.path.join(os.getcwd(), "DATA", "FIMO_OUT", "background_transcriptome_pval1e-2")

if pval_1000:
    diagram_title = "Enrichment of autologous binding via FIMO\nFull Transcriptome background of first order markov chain\nP-value cutoff: 1e-3"
    data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "background_transcriptome_pval1e-3")

if pval_10000:
    if transcriptome_background:
        diagram_title = "Enrichment of autologous binding via FIMO\nFull Transcriptome background of first order markov chain\nP-value cutoff: 1e-3"
        data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "background_transcriptome")

    if individual_background:
        diagram_title = "Enrichment of autologous binding via FIMO\nIndividual transcript parts as background - 1st order markov chain\nP-value cutoff: 1e+4"
        data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "background_individual")


directories = os.listdir(data_path)  # list of directiories containing output tsv files (& other files)


experiments = ["RNAcompete", "SELEX", "HT-SELEX"]
subsequences = ["UTR3", "UTR5", "CDS", "transcript"]


rnacompete_ppms = list(attract_ppms["RNAcomp"].keys( ))
selex_ppms = list(attract_ppms["SELEX"].keys( ))
htselex_ppms = list(htselex_ppms.keys( ))
rbps = [rnacompete_ppms, selex_ppms, htselex_ppms]


def rename_MANE_seqs_duplicated_motifs(MANE_transcriptome, rnacompete_ppms):
    for ppm in rnacompete_ppms:
        if ppm not in MANE_transcriptome:
            ppm_name_without_extension = str(ppm)[:-2]
            MANE_transcriptome[ppm] = MANE_transcriptome[ppm_name_without_extension]

rename_MANE_seqs_duplicated_motifs(MANE_transcriptome, rnacompete_ppms)

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
        header = f.readline()
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

    for key in MANE_transcriptome.keys( ):
        MANE_transcriptome[key]["transcript"] = MANE_transcriptome[key].pop("cDNA")

    # great outer loop for going through files
    for i, exp in enumerate(experiments):
        for subseq in subsequences:


            tsv_file_path = matches_sorted[exp][subseq]
            infos = read_tsv_file(tsv_file_path)


            add_sequence_length_to_infos(infos, subseq)
            print(f">>> SEQUENCE LENGTHS HAVE BEEN ADDED TO THE DICTIONARY\n")


            group_by_motif_id_and_sequence_id(infos)


            calc_coverages_autol_len(infos,
                                     subseq,
                                     len_normalize=False,
                                     consider_overlap=False,
                                     background_longer_than_autol_only=False)
            print(f">>> DONE WITH COVERAGE CALCULATIONS\n")

            autologous, background = []
            for motif in infos.keys():

                mean_cov_per_motif, std_cov_per_motif = get_mean_std(infos[motif])

                autologous, background = calc_z_scores(infos[motif],
                                                    mean_cov_per_motif,
                                                    std_cov_per_motif,
                                                    motif)

                autologous.append(autologous)
                background.append(background)

            print(f">>> Z-SCORES HAVE BEEN CALCULATED\n")


            p_val = binary_vs_averaged(autologous, background, ranked="no", side="higher")

            number_of_motifs_used = count_amount_of_motifs_per_experiment(attract_ppms, htselex_ppms)

            plot_analysis_results(autologous,
                                  background,
                                  p_val,
                                  number_of_motifs_used,
                                  i, # see for-loop above: i helps put labels on the plot
                                  calc_transcript = True,
                                  calc_CDS = True,
                                  calc_UTR3 = True,
                                  calc_UTR5 = True)



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



def group_by_motif_id_and_sequence_id(matches_by_motif):
    motif_subseq = {}
    for match in matches_by_motif:
        seq_id = match[0]
        motif_id = match[1]
        if motif_id not in motif_subseq:
            motif_subseq[motif_id] = {}

        if seq_id not in motif_subseq[motif_id]:
            motif_subseq[motif_id][seq_id] = []

        motif_subseq[motif_id][seq_id].append(match)
    return motif_subseq




def calc_coverages_autol_len(matches_by_subseq, subseq, len_normalize=False, consider_overlap=False, background_longer_than_autol_only=False):

    for motif in matches_by_subseq:

        for seq in matches_by_subseq[motif].values(): #all the matches a motif has with a sequence

            autologous_seq_len = len(MANE_transcriptome[motif][subseq]) # get the length of the autologous sequence
            seq_len = seq[0][4]
            seq_id = seq[0][0]
            motif_id = seq[0][1]

            if background_longer_than_autol_only:

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

                if len_normalize:
                    cov[2] = cov[2] / seq_len
                    matches_by_subseq[motif].append(cov)
                else:
                    matches_by_subseq[motif].append(cov)



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

                if len_normalize:
                    cov[2] = cov[2] / seq_len
                    matches_by_subseq[motif].append(cov)
                else:
                    matches_by_subseq[motif].append(cov)



    #ppms[exp] = len(matches_by_subseq.keys())

    return matches_by_subseq, len_normalize, consider_overlap
# coverages = data; len_norm and consider_overlap are for specifying the mode of analysis during plotting;
# ppms gives the number of matching motifs to the plot (N = )
    # coverages: holds motifs and the sequences each motif matched with. For every sequence, there's a coverage



def get_mean_std(matches_by_sequences):

    mean_cov = []
    std_cov = []

    for box in matches_by_sequences( ):
        #box = all sequences a motif matched with
        for i in range(len(box)):
            mean_cov.append(box[i][2])
            std_cov.append(box[i][2])

    mean_cov = np.nanmean(mean_cov[i])
    std_cov = np.nanstd(std_cov)

    return mean_cov, std_cov


def build_up_final_distribution_dict(exp, subseq):
    sorter_subsequences = {"UTR3": "autologous UTR3",
                     "UTR5": "autologous UTR5",
                     "CDS": "autologous CDS",
                     "transcript": "autologous transcript"}
    subseq_alt = sorter_subsequences[subseq]
    final_distributions = {}
    final_distributions[exp] = {}
    final_distributions[exp][subseq] = []
    final_distributions[exp][subseq_alt] = []

    return final_distributions


def calc_z_scores(coverages_by_motif, mean_cov, std_cov, motif_id):

    mean = mean_cov
    std = std_cov
    background = []
    autologous = []

    for box in coverages_by_motif:
        autologous_match_occurred = False # variable to include motifs that had matches (at all) but not the autologous one

        # a "box" is a container. It's a list of lists containing all sequence matches
        # of a single motif.
        # box[0] would show all the matches the given motif (motif is the key in this nesting of the dict)
        # has with a single sequence. The matches themselves are lists of the form
        # [sequence id, motif id, coverage]

        # when it did have a match, then it contains at least 1 list of shape [sequence id, motif id, coverage],
        # where motif id is always the same and sequence ids are all sequences the motif had matches with

        for i in range(len(box)):
            seq_id = box[i][0]
            z_val = (box[i][2] - mean) / std
            background.append(z_val)

            if seq_id == motif_id:
                autologous.append(z_val)
                autologous_match_occurred = True


        if not autologous_match_occurred:
            autologous.append((0 - mean) / std)


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
def plot_analysis_results(autologous, background, pvalue, ppms, i, calc_transcript = True, calc_CDS = True, calc_UTR3 = True, calc_UTR5 = True):
    print("\n>>> PLOTTING THE DATA\n")

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

    fig = plt.figure(figsize=[18.3/2.54, 11.0/2.54], constrained_layout=True, dpi=300)
    grid = fig.add_gridspec(1,1)
    ax3 = plt.subplot(grid[0,0])

    ##########################################################
    # PLOT LEFT SIDE of individual diagrams (autologous part):
    ##########################################################
    if i == 0:
        if calc_transcript:
            ax3.scatter([i-0.05]*len(autologous), autologous, linewidths=0, color="red", alpha=0.2, zorder=3, label="autol. RNA (z-score in resp. transcriptome)")
            ax3.plot([i-0.25,i-0.05],[np.mean(autologous)]*2, color="red", zorder=2, label="(autologous) RNA (mean)", alpha=0.4)

        if calc_CDS:
            ax3.scatter([i-0.10]*len(autologous), autologous, linewidths=0, color="orange", alpha=0.2, zorder=3, label="autol. CDS (z-score in resp. transcriptome CDS)")
            ax3.plot([i-0.25,i-0.05],[np.mean(autologous)]*2, color="orange", zorder=2, label="(autologous) CDS (mean)", alpha=0.4)

        if calc_UTR3:
            ax3.scatter([i-0.15]*len(autologous), autologous, linewidths=0, color="lightblue", alpha=0.2, zorder=3, label="autol. UTR3 (z-score in resp. transcriptome UTR3)")
            ax3.plot([i-0.25,i-0.05],[np.nanmean(autologous)]*2, color="lightblue", zorder=2, label="(autologous) UTR3 (mean)", alpha=0.4)

        if calc_UTR5:
            ax3.scatter([i-0.20]*len(autologous), autologous, linewidths=0, color="blue", alpha=0.2, zorder=3, label="autol. UTR5 (z-score in resp. transcriptome UTR5)")
            ax3.plot([i-0.25,i-0.05],[np.nanmean(autologous)]*2, color="blue", zorder=2, label="(autologous) UTR5 (mean)", alpha=0.4)

    else:
        if calc_transcript:
            ax3.scatter([i-0.05]*len(autologous), autologous, linewidths=0, color="red", alpha=0.2, zorder=3)
            ax3.plot([i-0.25,i-0.05],[np.mean(autologous)]*2, color="red", zorder=2, alpha=0.4)
        if calc_CDS:
            ax3.scatter([i-0.10]*len(autologous), autologous, linewidths=0, color="orange", alpha=0.2, zorder=3)
            ax3.plot([i-0.25,i-0.05],[np.mean(autologous)]*2, color="orange", zorder=2, alpha=0.4)
        if calc_UTR3:
            ax3.scatter([i-0.15]*len(autologous), autologous, linewidths=0, color="lightblue", alpha=0.2, zorder=3)
            ax3.plot([i-0.25,i-0.05],[np.mean(autologous)]*2, color="lightblue", zorder=2, alpha=0.4)

        if calc_UTR5:
            ax3.scatter([i-0.20]*len(autologous), autologous, linewidths=0, color="blue", alpha=0.2, zorder=3)
            ax3.plot([i-0.25,i-0.05],[np.nanmean(autologous)]*2, color="blue", zorder=2, alpha=0.4)


    if calc_transcript:
        hist = vertical_hist(background, i+0.0, 0.6, borders=[-6,6], bin_number="standard", window_average=5)
        vxs_trans, vys_trans, mean_trans = hist["x_smooth"],hist["y_smooth"],hist["mean"]

    if calc_CDS:
        hist = vertical_hist(background, i + 0.0, 0.6, borders=[-6, 6], bin_number="standard", window_average=5)
        vxs_CDS, vys_CDS, mean_CDS = hist["x_smooth"],hist["y_smooth"],hist["mean"]

    if calc_UTR3:
        hist = vertical_hist(background, i+0.0, 0.6, borders=[-6,6], bin_number="standard", window_average=5)
        vxs_UTR3, vys_UTR3, mean_UTR3 = hist["x_smooth"],hist["y_smooth"],hist["mean"]

    if calc_UTR5:
        hist = vertical_hist(background, i+0.0, 0.6, borders=[-6,6], bin_number="standard", window_average=5)
        vxs_UTR5, vys_UTR5, mean_UTR5 = hist["x_smooth"],hist["y_smooth"],hist["mean"]


    ##############################################################
    # PLOT RIGHT SIDE of individual diagrams (transcriptome part):
    ##############################################################
    if i == 0:
        if calc_transcript:
            ax3.plot(vxs_trans,vys_trans, color="red", zorder=2, label="transcriptome(s)\n(merged histogram)", alpha=0.5, linewidth=1)
            ax3.plot([i+0.05,i+0.25],[mean_trans]*2, color="red", zorder=2, alpha=0.5)
        if calc_CDS:
            ax3.plot(vxs_CDS,vys_CDS, color="orange", zorder=2, label="CDS (merged histogram)", alpha=0.5, linewidth=1)
            ax3.plot([i+0.05,i+0.25],[mean_CDS]*2, color="orange", zorder=2, alpha=0.5)
        if calc_UTR3:
            ax3.plot(vxs_UTR3,vys_UTR3, color="lightblue", zorder=2, label="UTR3 (merged histogram)", alpha=0.5, linewidth=1)
            ax3.plot([i+0.05,i+0.25],[mean_UTR3]*2, color="lightblue", zorder=2, alpha=0.5)
        if calc_UTR5:
            ax3.plot(vxs_UTR5,vys_UTR5, color="blue", zorder=2, label="UTR5 (merged histogram)", alpha=0.5, linewidth=1)
            ax3.plot([i+0.05,i+0.25],[mean_UTR5]*2, color="blue", zorder=2, alpha=0.5)

    else:
        if calc_transcript:
            ax3.plot(vxs_trans,vys_trans, color="red", zorder=2, alpha=0.5, linewidth=1)
            ax3.plot([i+0.05,i+0.25],[mean_trans]*2, color="red", zorder=2, alpha=0.5)
        if calc_CDS:
            ax3.plot(vxs_CDS,vys_CDS, color="orange", zorder=2, alpha=0.5, linewidth=1)
            ax3.plot([i+0.05,i+0.25],[mean_CDS]*2, color="orange", zorder=2, alpha=0.5)
        if calc_UTR3:
            ax3.plot(vxs_UTR3,vys_UTR3, color="lightblue", zorder=2, alpha=0.5, linewidth=1)
            ax3.plot([i+0.05,i+0.25],[mean_UTR3]*2, color="lightblue", zorder=2, alpha=0.5)
        if calc_UTR5:
            ax3.plot(vxs_UTR5,vys_UTR5, color="blue", zorder=2, alpha=0.5, linewidth=1)
            ax3.plot([i+0.05,i+0.25],[mean_UTR5]*2, color="blue", zorder=2, alpha=0.5)

    # Write number of RBPs and p values for different experiments:
    ax3.text(i, -6.5, "N = %.0f"%ppms[i], ha="center")
    if calc_transcript:
        ax3.text(i, -7, "$p_{transcript}$ = %.7f"%pvalue, ha="center", color="red")
    if calc_CDS:
        ax3.text(i, -7.5, "$p_{CDS}$ = %.7f"%pvalue, ha="center", color="orange")
    if calc_UTR3:
        ax3.text(i, -8, "$p_{UTR3}$ = %.7f"%pvalue, ha="center", color="lightblue")
    if calc_UTR5:
        ax3.text(i, -8.5, "$p_{UTR5}$ = %.7f"%pvalue, ha="center", color="blue")

    ax3.set_ylim([-6,6])
    ticklist = [i for i in range(len(my_EXPERIMENTS))]
    ax3.xaxis.tick_top()
    ax3.set_xticks(ticklist)
    ax3.set_xticklabels(my_EXPERIMENTS)

    ax3.legend(ncol=2, edgecolor="grey", framealpha = 0.8, facecolor = "white", frameon=True, loc=4)
    ax3.set_ylabel("motif coverage (nt/nt), z-score")

    global diagram_title
    #diagram_title = "FIMO analysis with p-value cutoff of 1e-4 - coverages length normalized"

    ax3.set_title(diagram_title)

    plt.grid()
    plt.show()

pipeline_for_FIMO_analysis(matches_sorted)