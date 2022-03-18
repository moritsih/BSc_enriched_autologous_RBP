def fetch_seq_len(matches_sorted_dict):
    global MANE_transcriptome
    global experiments
    global subsequences


    print(">>> MANE DICT KEY HAS BEEN CHANGED\n")

    for exp in experiments:
        for subseq in subsequences:
            for i in range(len(matches_sorted_dict[exp][subseq])):
                seq_id = matches_sorted_dict[exp][subseq][i][0]
                seq_len = len(MANE_transcriptome[seq_id][subseq])
                matches_sorted_dict[exp][subseq][i].append(seq_len)



    return matches_sorted_dict


def get_set_of_motifs(matches):
    global experiments
    global subsequences

# Nested dictionary is initiated with all necessary keys
    motifs = {}

    for exp in matches.keys():
        motif_exp = {}
        motifs[exp] = motif_exp

        for subseq in subsequences:
            motif_exp[subseq] = group_by_motif_id_and_sequence_id(matches[exp][subseq])

        print(f">>> {exp} MOTIFS ARE NOW IN A NICE LOOK-UP TABLE\n")

    return motifs


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
