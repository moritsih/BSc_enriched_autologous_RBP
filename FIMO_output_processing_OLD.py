import os
from matplotlib import pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy import stats
import random
from FIMO_input_processing import MANE_transcriptome, attract_ppms, htselex_ppms, dict_MANE
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



def fetch_matches():
    global data_path
    global directories

    matches_files = {}

    for dir in directories:  # loop through folders of different experiment/subsequence combinations
        matches_files[dir] = []
        exp = os.path.join(data_path, dir)  # define which experiment we're dealing with
        for file in os.listdir(exp):
            if ".tsv" in file:  # there are other types of files in the dir
                with open(os.path.join(exp, file), "r") as f:
                    _ = f.readline( )
                    content = f.read( ).split("\n")
                    content = [x for x in content if x and not x.startswith("#")]  # removing bottom lines
                    for line in content:
                        line = line.split("\t")
                        motif_id = line[0]
                        seq_id = line[2]
                        start = line[3]
                        stop = line[4]
                        info = [seq_id, motif_id, start, stop]
                        matches_files[dir].append(info)

        print(f">>> FETCHED MATCHES FROM {dir}\n")

    return matches_files


matches_raw = fetch_matches( )


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
matches_raw = None

def plot_num_matches(match_dict):
    global individual_background, transcriptome_background

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    if individual_background:
        fig.suptitle("Number of matches - background distributions by subsequences")
    if transcriptome_background:
        fig.suptitle("Number of matches - transcriptome background distribution")

    for i, exp in enumerate(matches_sorted.keys( )):
        subseq_matches = []
        seqs = list(matches_sorted[exp].keys( ))
        for seq in seqs:
            lengths = len(matches_sorted[exp][seq])
            subseq_matches.append(lengths)
        sns.set_style=("dark")
        sns.barplot(ax=axes[i], x=seqs, y=subseq_matches)
        axes[i].set_title(exp)

    plt.show( )

#plot_num_matches(matches_sorted)


def fetch_seq_len(matches_sorted_dict):
    global MANE_transcriptome
    global experiments
    global subsequences

    for key in MANE_transcriptome.keys( ):
        dict_MANE[key]["transcript"] = dict_MANE[key].pop("cDNA")
        MANE_transcriptome[key]["transcript"] = MANE_transcriptome[key].pop("cDNA")

    print(">>> MANE DICT KEY HAS BEEN CHANGED\n")

    for exp in experiments:
        for subseq in subsequences:
            for i in range(len(matches_sorted_dict[exp][subseq])):
                seq_id = matches_sorted_dict[exp][subseq][i][0]
                seq_len = len(MANE_transcriptome[seq_id][subseq])
                matches_sorted_dict[exp][subseq][i].append(seq_len)

        print(f">>> SEQUENCE LENGTHS OF {exp} HAVE BEEN ADDED TO THE DICTIONARY\n")

    return matches_sorted_dict


matches = fetch_seq_len(matches_sorted)
matches_sorted = None

def get_set_of_motifs(matches):
    global experiments
    global subsequences

# Nested dictionary is initiated with all necessary keys
    motifs = {}
    for exp in experiments:
        motifs[exp] = {}
        for subseq in subsequences:
            motifs[exp][subseq] = {}

    for exp in matches.keys():
        for subseq in subsequences:
            motif_set = [] # list holding all motifs that had matches (usually less than full scope of motifs)
            for match in matches[exp][subseq]:
                motif_set.append(match[1]) # motif ID is appended
            motif_set = set(motif_set) # need each matching motif only once, set() operation
            for motif in motif_set:
                motifs[exp][subseq][motif] = {} # initiate sub-dictionary for each matching motif


    for exp in matches.keys( ):
        for subseq in matches[exp].keys( ):
            for match in matches[exp][subseq]: # all matches that happened
                motif_id = match[1] #
                seq_id = match[0]
                try:
                    motifs[exp][subseq][motif_id][seq_id].append(match)
                except:
                    motifs[exp][subseq][motif_id][seq_id] = []
                    motifs[exp][subseq][motif_id][seq_id].append(match)

        print(f">>> {exp} MOTIFS ARE NOW IN A NICE LOOK-UP TABLE\n")

    return motifs

matches = get_set_of_motifs(matches)


def calc_coverages_autol_len(len_normalize=False, consider_overlap=False):
    global matches, experiments

    coverages = {}
    ppms = {}
    for exp in matches.keys( ):
        coverages[exp] = {}
        for subseq in matches[exp].keys( ):
            coverages[exp][subseq] = {}

            for motif in matches[exp][subseq].keys( ):
                coverages[exp][subseq][motif] = []

                for seq in matches[exp][subseq][motif].values(): #all the matches a motif has with a sequence

                    autol_len = len(dict_MANE[motif][subseq]) # get the length of the autologous sequence
                    seq_len = seq[0][4]

                    if seq_len >= autol_len:

                        for idx in range(len(seq)): # how many matches are there with a single sequence? idx will tell you

                            data = seq[idx]
                            seq_id = data[0]
                            motif_id = data[1]
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

                        start_idx_range = seq_len-autol_len
                        if start_idx_range == 0:
                            cov = [seq_id, motif_id, np.mean(array)]
                            continue
                        start_idx = np.random.randint(0, start_idx_range)
                        stop_idx = start_idx + autol_len

                        array = array[start_idx:stop_idx]

                        cov = [seq_id, motif_id, np.mean(array)]

                    else:
                        continue

                    if len_normalize:
                        cov[2] = cov[2] / seq_len
                        coverages[exp][subseq][motif].append(cov)
                    else:
                        coverages[exp][subseq][motif].append(cov)

        print(f">>> DONE WITH COVERAGE CALCULATIONS OF {exp}\n")

    for exp in experiments:
        for subseq in subsequences:
            ppms[exp] = len(coverages[exp][subseq].keys())

    return coverages, len_normalize, consider_overlap, ppms
    # coverages = data; len_norm and consider_overlap are for specifying the mode of analysis during plotting;
    # ppms gives the number of matching motifs to the plot (N = )
    # coverages: holds motifs and the sequences each motif matched with. For every sequence, there's a coverage

coverages, len_normalized, overlaps_considered, ppms = calc_coverages_autol_len(len_normalize=False, consider_overlap=False)


########################################################################################################################
def calc_coverages(matches_sorted_dict, len_normalize=False, consider_overlap=False):
    global matches, experiments

    coverages = {}
    ppms = {}
    for exp in matches_sorted_dict.keys():
        coverages[exp] = {}
        for subseq in matches_sorted_dict[exp].keys():
            coverages[exp][subseq] = {}

            for motif in matches_sorted_dict[exp][subseq].keys():
                coverages[exp][subseq][motif] = []

                for seq in matches_sorted_dict[exp][subseq][motif].values():
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
                        coverages[exp][subseq][motif].append(cov)
                    else:
                        coverages[exp][subseq][motif].append(cov)

        print(f">>> DONE WITH COVERAGE CALCULATIONS OF {exp}\n")

    for exp in experiments:
        for subseq in subsequences:
            ppms[exp] = len(coverages[exp][subseq].keys())

    return coverages, len_normalize, consider_overlap, ppms

#coverages, len_normalized, overlaps_considered, ppms = calc_coverages(matches, len_normalize=False, consider_overlap=False)
matches = None
#######################################################################################################################



def get_mean_std(coverages_dict):
    global experiments, subsequences

    mean_cov = {}
    std_cov = {}

    for exp in experiments:
        mean_cov[exp] = {}
        std_cov[exp] = {}

        for subseq in subsequences:
            mean_cov[exp][subseq] = {}
            std_cov[exp][subseq] = {}

            for motif_id, box in coverages_dict[exp][subseq].items( ):
                mean_cov[exp][subseq][motif_id] = []
                std_cov[exp][subseq][motif_id] = []
                #box = all sequences a motif matched with
                for i in range(len(box)):
                    mean_cov[exp][subseq][motif_id].append(box[i][2])
                    std_cov[exp][subseq][motif_id].append(box[i][2])

    for exp in experiments:
        for subseq in subsequences:
            for motif in mean_cov[exp][subseq].keys():
                mean_cov[exp][subseq][motif] = np.nanmean(mean_cov[exp][subseq][motif])
                std_cov[exp][subseq][motif] = np.nanstd(std_cov[exp][subseq][motif])

    return mean_cov, std_cov

mean_cov, std_cov = get_mean_std(coverages)


def calc_z_scores(dict_coverages, mean_cov, std_cov):
    global experiments, subsequences

    # dictionary used for creating the keys of autologous matches
    sorter_subseq = {"UTR3": "autologous UTR3",
                     "UTR5": "autologous UTR5",
                     "CDS": "autologous CDS",
                     "transcript": "autologous transcript"}

    final_distributions = {}

    for exp in experiments:
        final_distributions[exp] = {}

        for subseq in subsequences:

            subseq_a = sorter_subseq[subseq]
            final_distributions[exp][subseq_a] = []
            final_distributions[exp][subseq] = []



            for motif, box in dict_coverages[exp][subseq].items( ):
                mean = mean_cov[exp][subseq][motif]
                std = std_cov[exp][subseq][motif]
                autologous_match_occurred = False # variable to include motifs that had matches (at all) but not the autologous one
                background = []
                # a "box" is a container. It's a list of lists containing all sequence matches
                # of a single motif.
                # box[0] would show all the matches the given motif (motif is the key in this nesting of the dict)
                # has with a single sequence. The matches themselves are lists of the form
                # [sequence id, motif id, coverage]

                # when it did have a match, then it contains at least 1 list of shape [sequence id, motif id, coverage],
                # where motif id is always the same and sequence ids are all sequences the motif had matches with

                for i in range(len(box)):
                    seq_id = box[i][0]
                    motif_id = box[i][1]
                    z_val = (box[i][2] - mean) / std
                    background.append(z_val)
            # in order to get the right "shape" of data for the p-value function, I need
            # final_distributions[exp]["transcript"] = [[background], [background],...]

                    if seq_id == motif_id:
                        final_distributions[exp][subseq_a].append(z_val)

                        autologous_match_occurred = True


                if not autologous_match_occurred:
                    final_distributions[exp][subseq_a].append((0 - mean) / std)

                final_distributions[exp][subseq].append(background)


        print(f">>> Z-SCORES FOR {exp} HAVE BEEN CALCULATED\n")

    return final_distributions

final_distributions = calc_z_scores(coverages, mean_cov, std_cov)
coverages = None


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
def plot_analysis_results():
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

    calc_transcript = True
    calc_CDS = True
    calc_UTR3 = True
    calc_UTR5 = True

    for i, experiment in enumerate(my_EXPERIMENTS):
        ##########################################################
        # PLOT LEFT SIDE of individual diagrams (autologous part):
        ##########################################################
        if i == 0:
            if calc_transcript:
                #print(len(final_distributions[experiment]["transcript"]))
                #pr_trans_mean = (np.mean(final_distributions[experiment]["transcript"]))*2
                #pr_trans_mean_auto = (np.mean(final_distributions[experiment]["autologous transcript"]))*2
                #print(f"transcriptome mean:     {pr_trans_mean}")
                #print(f"autol. transcript mean: {pr_trans_mean_auto}")
                #outputfile.write(f"\ntranscriptome mean: {pr_trans_mean}")
                #outputfile.write(f"\nautol. transcript mean: {pr_trans_mean_auto}")

                ax3.scatter([i-0.05]*len(final_distributions[experiment]["autologous transcript"]), final_distributions[experiment]["autologous transcript"], linewidths=0, color="red", alpha=0.2, zorder=3, label="autol. RNA (z-score in resp. transcriptome)")
                ax3.plot([i-0.25,i-0.05],[np.mean(final_distributions[experiment]["autologous transcript"])]*2, color="red", zorder=2, label="(autologous) RNA (mean)", alpha=0.4)
            if calc_CDS:

                #pr_CDS_mean = (np.mean(final_distributions[experiment]["CDS"]))*2
                #pr_CDS_mean_auto = (np.mean(final_distributions[experiment]["autologous CDS"]))*2
                #print(f"CDS mean:     {pr_CDS_mean}")
                #print(f"autol. CDS mean: {pr_CDS_mean_auto}")
                #outputfile.write(f"\nCDS mean: {pr_CDS_mean}")
                #outputfile.write(f"\nautol. CDS mean: {pr_CDS_mean_auto}")

                ax3.scatter([i-0.10]*len(final_distributions[experiment]["autologous CDS"]), final_distributions[experiment]["autologous CDS"], linewidths=0, color="orange", alpha=0.2, zorder=3, label="autol. CDS (z-score in resp. transcriptome CDS)")
                ax3.plot([i-0.25,i-0.05],[np.mean(final_distributions[experiment]["autologous CDS"])]*2, color="orange", zorder=2, label="(autologous) CDS (mean)", alpha=0.4)
            if calc_UTR3:

                #pr_UTR3_mean = (np.nanmean(final_distributions[experiment]["UTR3"]))*2
                #pr_UTR3_mean_auto = (np.nanmean(final_distributions[experiment]["autologous UTR3"]))*2
                #print(f"UTR3 mean:     {pr_UTR3_mean}")
                #print(f"autol. UTR3 mean: {pr_UTR3_mean_auto}")
                #outputfile.write(f"\nUTR3 mean: {pr_UTR3_mean}")
                #outputfile.write(f"\nautol. UTR3 mean: {pr_UTR3_mean_auto}")

                ax3.scatter([i-0.15]*len(final_distributions[experiment]["autologous UTR3"]), final_distributions[experiment]["autologous UTR3"], linewidths=0, color="lightblue", alpha=0.2, zorder=3, label="autol. UTR3 (z-score in resp. transcriptome UTR3)")
                ax3.plot([i-0.25,i-0.05],[np.nanmean(final_distributions[experiment]["autologous UTR3"])]*2, color="lightblue", zorder=2, label="(autologous) UTR3 (mean)", alpha=0.4)
            if calc_UTR5:

                #pr_UTR5_mean = (np.nanmean(final_distributions[experiment]["UTR5"]))*2
                #pr_UTR5_mean_auto = (np.nanmean(final_distributions[experiment]["autologous UTR5"]))*2
                #print(f"UTR5 mean:     {pr_UTR5_mean}")
                #print(f"autol. UTR5 mean: {pr_UTR5_mean_auto}")
                #outputfile.write(f"\nUTR5 mean: {pr_UTR5_mean}")
                #outputfile.write(f"\nautol. UTR5 mean: {pr_UTR5_mean_auto}")

                ax3.scatter([i-0.20]*len(final_distributions[experiment]["autologous UTR5"]), final_distributions[experiment]["autologous UTR5"], linewidths=0, color="blue", alpha=0.2, zorder=3, label="autol. UTR5 (z-score in resp. transcriptome UTR5)")
                ax3.plot([i-0.25,i-0.05],[np.nanmean(final_distributions[experiment]["autologous UTR5"])]*2, color="blue", zorder=2, label="(autologous) UTR5 (mean)", alpha=0.4)

        else:
            if calc_transcript:
                #pr_trans_mean = (np.mean(final_distributions[experiment]["transcript"]))#*2
                #pr_trans_mean_auto = (np.mean(final_distributions[experiment]["autologous transcript"]))#*2
                #print(f"transcriptome mean:     {pr_trans_mean}")
                #print(f"autol. transcript mean: {pr_trans_mean_auto}")
                #outputfile.write(f"\ntranscriptome mean: {pr_trans_mean}")
                #outputfile.write(f"\nautol. transcript mean: {pr_trans_mean_auto}")

                ax3.scatter([i-0.05]*len(final_distributions[experiment]["autologous transcript"]), final_distributions[experiment]["autologous transcript"], linewidths=0, color="red", alpha=0.2, zorder=3)
                ax3.plot([i-0.25,i-0.05],[np.mean(final_distributions[experiment]["autologous transcript"])]*2, color="red", zorder=2, alpha=0.4)
            if calc_CDS:
                #pr_CDS_mean = (np.mean(final_distributions[experiment]["CDS"]))*2
                #pr_CDS_mean_auto = (np.mean(final_distributions[experiment]["autologous CDS"]))*2
                #print(f"CDS mean:     {pr_CDS_mean}")
                #print(f"autol. CDS mean: {pr_CDS_mean_auto}")
                #outputfile.write(f"\nCDS mean: {pr_CDS_mean}")
                #outputfile.write(f"\nautol. CDS mean: {pr_CDS_mean_auto}")

                ax3.scatter([i-0.10]*len(final_distributions[experiment]["autologous CDS"]), final_distributions[experiment]["autologous CDS"], linewidths=0, color="orange", alpha=0.2, zorder=3)
                ax3.plot([i-0.25,i-0.05],[np.mean(final_distributions[experiment]["autologous CDS"])]*2, color="orange", zorder=2, alpha=0.4)
            if calc_UTR3:
                #pr_UTR3_mean = (np.mean(final_distributions[experiment]["UTR3"]))*2
                #pr_UTR3_mean_auto = (np.mean(final_distributions[experiment]["autologous UTR3"]))*2
                #print(f"UTR3 mean:     {pr_UTR3_mean}")
                #print(f"autol. UTR3 mean: {pr_UTR3_mean_auto}")
                #outputfile.write(f"\nUTR3 mean: {pr_UTR3_mean}")
                #outputfile.write(f"\nautol. UTR3 mean: {pr_UTR3_mean_auto}")

                ax3.scatter([i-0.15]*len(final_distributions[experiment]["autologous UTR3"]), final_distributions[experiment]["autologous UTR3"], linewidths=0, color="lightblue", alpha=0.2, zorder=3)
                ax3.plot([i-0.25,i-0.05],[np.mean(final_distributions[experiment]["autologous UTR3"])]*2, color="lightblue", zorder=2, alpha=0.4)
            if calc_UTR5:
                #pr_UTR5_mean = (np.nanmean(final_distributions[experiment]["UTR5"]))*2
                #pr_UTR5_mean_auto = (np.nanmean(final_distributions[experiment]["autologous UTR5"]))*2
                #print(f"UTR5 mean:     {pr_UTR5_mean}")
                #print(f"autol. UTR5 mean: {pr_UTR5_mean_auto}")
                #outputfile.write(f"\nUTR5 mean: {pr_UTR5_mean}")
                #outputfile.write(f"\nautol. UTR5 mean: {pr_UTR5_mean_auto}")

                ax3.scatter([i-0.20]*len(final_distributions[experiment]["autologous UTR5"]), final_distributions[experiment]["autologous UTR5"], linewidths=0, color="blue", alpha=0.2, zorder=3)
                ax3.plot([i-0.25,i-0.05],[np.nanmean(final_distributions[experiment]["autologous UTR5"])]*2, color="blue", zorder=2, alpha=0.4)


        print(f">>> RANDOMIZING {experiment} DATA\n")
        # TRANSCRIPTOME DATA
        if calc_transcript:
            p_trans = binary_vs_averaged(final_distributions[experiment]["autologous transcript"],final_distributions[experiment]["transcript"],ranked="no")

            hist = vertical_hist(final_distributions[experiment]["transcript"], i+0.0, 0.6, borders=[-6,6], bin_number="standard", window_average=5)
            vxs_trans, vys_trans, mean_trans = hist["x_smooth"],hist["y_smooth"],hist["mean"]


            #print("mean_trans", mean_trans)
            #outputfile.write(f"\nmean_trans: {mean_trans}")
        # CDS DATA
        if calc_CDS:
            p_CDS = binary_vs_averaged(final_distributions[experiment]["autologous CDS"],final_distributions[experiment]["CDS"],ranked="no")
            hist = vertical_hist(final_distributions[experiment]["CDS"], i+0.0, 0.6, borders=[-6,6], bin_number="standard", window_average=5)
            vxs_CDS, vys_CDS, mean_CDS = hist["x_smooth"],hist["y_smooth"],hist["mean"]

            #print("mean_CDS", mean_CDS)
            #outputfile.write(f"\nmean_CDS: {mean_CDS}")
        # UTR3 DATA
        if calc_UTR3:
            p_UTR3 = binary_vs_averaged(final_distributions[experiment]["autologous UTR3"],final_distributions[experiment]["UTR3"],ranked="no")
            hist = vertical_hist(final_distributions[experiment]["UTR3"], i+0.0, 0.6, borders=[-6,6], bin_number="standard", window_average=5)
            vxs_UTR3, vys_UTR3, mean_UTR3 = hist["x_smooth"],hist["y_smooth"],hist["mean"]

            #x = np.hstack(final_distributions[experiment]["UTR3"])
            #x = x[np.isfinite(x)]
            #vert_hist = np.histogram(x, bins=100)
            #mean_UTR3 =  np.mean(vert_hist[1])
            #vxs_UTR3, vys_UTR3 = vert_hist[0], vert_hist[1][:-1]

            #print("mean_UTR3", mean_UTR3)
            #outputfile.write(f"\nmean_UTR3: {mean_UTR3}")
        # UTR5 DATA
        if calc_UTR5:
            p_UTR5 = binary_vs_averaged(final_distributions[experiment]["autologous UTR5"],final_distributions[experiment]["UTR5"],ranked="no")
            hist = vertical_hist(final_distributions[experiment]["UTR5"], i+0.0, 0.6, borders=[-6,6], bin_number="standard", window_average=5)
            vxs_UTR5, vys_UTR5, mean_UTR5 = hist["x_smooth"],hist["y_smooth"],hist["mean"]


            #print("mean_UTR5", mean_UTR5)
            #outputfile.write(f"\nmean_UTR5: {mean_UTR5}")


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
        ax3.text(i, -6.5, "N = %.0f"%ppms[experiment], ha="center")
        if calc_transcript:
            ax3.text(i, -7, "$p_{transcript}$ = %.7f"%p_trans, ha="center", color="red")
        if calc_CDS:
            ax3.text(i, -7.5, "$p_{CDS}$ = %.7f"%p_CDS, ha="center", color="orange")
        if calc_UTR3:
            ax3.text(i, -8, "$p_{UTR3}$ = %.7f"%p_UTR3, ha="center", color="lightblue")
        if calc_UTR5:
            ax3.text(i, -8.5, "$p_{UTR5}$ = %.7f"%p_UTR5, ha="center", color="blue")

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


plot_analysis_results()
