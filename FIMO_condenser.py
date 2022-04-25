import os
import sys
from additional_code.read_MANE import dict_MANE
from FIMO_input_processing import attract_ppms, htselex_ppms

ppms = [attract_ppms["RNAcomp"], attract_ppms["SELEX"], htselex_ppms]

# this script requires 1 commandline argument; for ex. enter "python3 FIMO_condenser.py pval1e-4"
pval_cutoff = sys.argv[1]

if pval_cutoff == "5e-2":
    data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "pval5e-2")

if pval_cutoff == "1e-2":
    data_path = os.path.join(os.getcwd(), "DATA", "FIMO_OUT", "pval1e-2")

if pval_cutoff == "1e-3":
    data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "pval1e-3")

if pval_cutoff == "1e-4":
    data_path = os.path.join(os.getcwd( ), "DATA", "FIMO_OUT", "pval1e-4")


files = [x for x in os.listdir(data_path) if "full" in x]
target_location_file = os.path.join(os.getcwd(), "DATA", "condensed_files")

def create_condensed_file_from_raw_matches(dict_MANE, data_path, files, ppms, target_location_file):

    for exp in ppms:
        rename_MANE_seqs_duplicated_motifs(dict_MANE, exp)

    match_files = store_filenames_for_retrieval(data_path, files)

    for file in match_files:

        match_per_motif_n_sequence_generator = group_by_motif_id_and_sequence_id(file)

        filename = path_to_filename_converter(file)

        ppm = filename_to_ppm_identifier_converter(filename, attract_ppms, htselex_ppms)

        with open(os.path.join(target_location_file, filename), "a") as f:
            header = "Protein_gene-name\t" \
                     "matrix_ID\t" \
                     "matrix_length\t" \
                     "transcript_gene-name\t" \
                     "transcript_ID\t" \
                     "UTR5_length(nt)\t" \
                     "CDS_length(nt)\t" \
                     "UTR3_length(nt)\n"

            f.write(header)

            for matches_per_sequence in match_per_motif_n_sequence_generator:

                write_header_per_motif_seq_pair(dict_MANE, f, matches_per_sequence, ppm)

                for single_match_per_seq in matches_per_sequence:

                    start = single_match_per_seq[2]
                    pval = single_match_per_seq[3]

                    f.write(f"{start}\t{pval}\n")

                f.write("\n")



def write_header_per_motif_seq_pair(dict_MANE, f, matches_per_sequence, ppm):
    seq_id = matches_per_sequence[0][0]
    matrix_id = matches_per_sequence[0][1]

    protein_gene_name, \
    transcript_gene_name, \
    transcript_id, \
    utr5_len, \
    cds_len, \
    utr3_len = MANE_database_information_retriever(seq_id, matrix_id, dict_MANE)

    matrix_len = get_matrix_length(matrix_id, ppm)

    f.write(f"{protein_gene_name}\t"
            f"{matrix_id}\t"
            f"{matrix_len}\t"
            f"{transcript_gene_name}\t"
            f"{transcript_id}\t"
            f"{utr5_len}\t"
            f"{cds_len}\t"
            f"{utr3_len}\n")


# stores access to folders where the FIMO-output data lies; for retrieval inside loop
def store_filenames_for_retrieval(data_path, files):

    matches_files = []

    for file in files:
        tsv_file_path = os.path.join(data_path, file)
        matches_files.append(tsv_file_path)

        print(f">>> FETCHED MATCHES FROM {file}\n")

    return matches_files



def rename_MANE_seqs_duplicated_motifs(MANE_transcriptome, ppms):
    for ppm in ppms:
        if ppm not in MANE_transcriptome:
            ppm_name_without_extension = str(ppm)[:-2]
            MANE_transcriptome[ppm] = MANE_transcriptome[ppm_name_without_extension]

    return MANE_transcriptome



def read_tsv_file(tsv_file_path):

    with open(tsv_file_path, "r") as f:
        _ = f.readline()
        content = f.read().split("\n")
        content = [x for x in content if x and not x.startswith("#")]  # removing bottom lines

        for line in content:
            line = line.split("\t")
            matrix_id = line[0]
            seq_id = line[2]
            start = line[3]
            pval = line[7]

            infos = [seq_id, matrix_id, start, pval]

            yield infos


def group_by_motif_id_and_sequence_id(file_path):

    matches_sorted = {}

    info_generator = read_tsv_file(file_path)
    for infos in info_generator:

        seq_id = infos[0]
        motif_id = infos[1]

        if motif_id not in matches_sorted:
            matches_sorted[motif_id] = {}

        if seq_id not in matches_sorted[motif_id]:
            matches_sorted[motif_id][seq_id] = []

        matches_sorted[motif_id][seq_id].append(infos)

    for motif in matches_sorted.keys():
        for sequence in matches_sorted[motif].keys():
            yield matches_sorted[motif][sequence]



def MANE_database_information_retriever(seq_id, matrix_id, dict_MANE):

    protein_info = dict_MANE[matrix_id]
    transcript_info = dict_MANE[seq_id]

    protein_name = protein_info["mane_gene_name"]

    transcript_gene_name = transcript_info["mane_gene_name"]
    transcript_id = transcript_info["mane_transcript_ID"]
    cds_len = len(transcript_info["CDS"])
    utr5_len = len(transcript_info["UTR5"])
    utr3_len = len(transcript_info["UTR3"])

    return protein_name, transcript_gene_name, transcript_id, utr5_len, cds_len, utr3_len


def get_matrix_length(matrix_id, ppms):
    return len(ppms[matrix_id])


def path_to_filename_converter(filepath):

    if filepath.endswith("rnacompete_full.tsv"):
        exp = "RNAcompete.txt"
    elif filepath.endswith("selex_full.tsv"):
        exp = "SELEX.txt"
    elif filepath.endswith("htselex_full.tsv"):
        exp = "HT-SELEX.txt"
    else:
        print("Wtf happened")

    return exp



def filename_to_ppm_identifier_converter(filename, attract_ppms, htselex_ppms):

    converter_dict = {
        "RNAcompete.txt":attract_ppms["RNAcomp"],
        "SELEX.txt":attract_ppms["SELEX"],
        "HT-SELEX.txt":htselex_ppms
    }

    return converter_dict[filename]




create_condensed_file_from_raw_matches(dict_MANE, data_path, files, ppms, target_location_file)


