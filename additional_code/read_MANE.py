""" This module reads information from the MANE Select files and converts them
    to a nested dictionary which can be loaded by other programs. """

# %%
###################################################################################################
# IMPORTS #########################################################################################
###################################################################################################

import string
import os


home = os.getcwd()
data_path = (os.path.join(home,"DATA"))

###################################################################################################
# FUNCTIONS #######################################################################################
###################################################################################################

def select_nest(dict_2,keys,last,count=0):
    if last == 0:
        return dict_2
    if count == 0:
        count += 1
    if count == last:
        return dict_2[keys[count-1]]
    else:
        return select_nest(dict_2[keys[count-1]],keys,last,count+1)


def add_dict_branch(dict_1,keys,end_type):
    for count, key in enumerate(keys):
        if count < len(keys)-1:
            if key not in select_nest(dict_1,keys,count).keys():
                select_nest(dict_1,keys,count)[key] = {}
        else:
            if key not in select_nest(dict_1,keys,count).keys():
                select_nest(dict_1,keys,count)[key] = end_type


###################################################################################################
# CREATE MANE DICTIONARY ##########################################################################
###################################################################################################

# Filter clinical IDs from MANE file by creating a set of them:
mane_clinical_IDs = set()

print("\n>>> FILTERING MANE CLINICAL IDs")
with open(data_path + "/MANE/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf", "r") as clinicalfile:
    for counter, line in enumerate(clinicalfile):
        line = line.split()
        feature = line[2]
        if feature == "gene":
            continue
        transcript_ID = line[11].split('"')[1]
        tags = line[16:]
        if '"MANE_Select";' in tags:
            mane_status = True
        elif '"MANE_Plus_Clinical";' in tags:
            mane_status = False
            mane_clinical_IDs.add(transcript_ID)
        else:
            mane_status = "Unknown"
            print(mane_status)


# Read MANE file and create dictionary with IDs and sequence:
dict_MANE = {}

print("\n>>> READING MANE FNA FILE")
with open(data_path + "/MANE/MANE.GRCh38.v0.95.select_ensembl_rna.fna", "r") as manefile:
#with open(home+"/DATA/MANE/babyMANE.fna", "r") as manefile:
    for counter, line in enumerate(manefile):
        # Check if line is a FASTA header and split it:
        if line.startswith(">"):
            line = line.split()

            # Ignore if the transcript is in the set of clinical IDs:
            if line[0][1:] in mane_clinical_IDs:
                skip_line = True
                continue
            else:
                skip_line = False

            # Add GENE ID to the dictionary:
            mane_gene_ID = line[3]
            if "gene:" not in mane_gene_ID or "." not in mane_gene_ID:
                raise ValueError("gene: or . not in mane_gene_ID")
            else:
                mane_gene_ID = mane_gene_ID[5:20]
                if mane_gene_ID in dict_MANE.keys():
                    del dict_MANE[mane_gene_ID]
                    raise ValueError("gene ID not unique: ", mane_gene_ID)
                add_dict_branch(dict_MANE,[mane_gene_ID, "mane_gene_ID"],mane_gene_ID)
                add_dict_branch(dict_MANE,[mane_gene_ID,"cDNA"],"")

            # Add GENE NAME to the dictionary:
            mane_gene_name = line[6]
            if "gene_symbol:" not in mane_gene_name:
                raise ValueError("gene_symbol: not in mane_gene_name")
            else:
                mane_gene_name = mane_gene_name[12:]
                add_dict_branch(dict_MANE,[mane_gene_ID,"mane_gene_name"],mane_gene_name)

            # Add TRANSCRIPT ID to the dictionary:
            mane_transcript_ID = line[0]
            if ">" not in mane_transcript_ID or "." not in mane_transcript_ID :
                raise ValueError("> or . not in mane_transcript_ID")
            else:
                mane_transcript_ID = mane_transcript_ID[1:]
                add_dict_branch(dict_MANE,[mane_gene_ID,"mane_transcript_ID"],mane_transcript_ID)

            # Add GENOMIC LENGTH to the dictionary:
            mane_genomic_position = line[2].split(":")
            if len(mane_genomic_position) != 6:
                raise ValueError("genomic position aberrant")
            else:
                mane_genomic_position = [int(x) for x in mane_genomic_position[3:5]]
                add_dict_branch(dict_MANE,[mane_gene_ID,"mane_genomic_length"],mane_genomic_position[1]-mane_genomic_position[0]+1)

        # Check if line only contains allowed characters:
        elif len((set(string.printable)-set("ATCG"))&set(line[0][:-1])) != 0:
            raise ValueError("Sequence contains character outside ATCG")

        # Add lines below the header (= cDNA) to the dictionary, skip for clinical IDs:
        else:
            if skip_line:
                continue
            dict_MANE[mane_gene_ID]["cDNA"] += line.rstrip()


# Convert MANE cDNA to RNA by replacing "T" with "U":
for gene in dict_MANE.keys():
    sequence = dict_MANE[gene]["cDNA"].replace("T", "U")
    dict_MANE[gene]["cDNA"] = sequence


# Add protein sequences to the dictionary:
print("\n>>> READING MANE FAA FILE")
with open(data_path + "/MANE/MANE.GRCh38.v0.95.select_ensembl_protein.faa", "r") as manefile:
    for counter, line in enumerate(manefile):
        # Check if line is a FASTA header and split it:
        if line.startswith(">"):
            line = line.split(" ")

            # Ignore if the transcript is in the set of clinical IDs:
            if line[4].split(":")[1] in mane_clinical_IDs:
                skip_line = True
                continue
            else:
                skip_line = False

            # Add GENE ID to the dictionary:
            mane_gene_ID = line[3]
            if "gene:" not in mane_gene_ID or "." not in mane_gene_ID:
                raise ValueError("gene: or . not in mane_gene_ID")
            else:
                mane_gene_ID = mane_gene_ID[5:20]
                if mane_gene_ID not in dict_MANE.keys():
                    raise ValueError("gene for protein missing in transcript dataset: ", mane_gene_ID)
                else:
                    if "mane_protein_ID" in dict_MANE[mane_gene_ID].keys():
                        del dict_MANE[mane_gene_ID]["protein"]
                        print("gene not unique: ", mane_gene_ID)
                add_dict_branch(dict_MANE,[mane_gene_ID,"protein"],"")

            # Compare GENE NAME with the dictionary:
            mane_gene_name = line[7]
            if "gene_symbol:" not in mane_gene_name:
                raise ValueError("gene_symbol: not in mane_gene_name")
            else:
                mane_gene_name = mane_gene_name[12:]
                if mane_gene_name != dict_MANE[mane_gene_ID]["mane_gene_name"]:
                    print(mane_gene_ID, mane_gene_name, dict_MANE[mane_gene_ID]["mane_gene_ID"])
                    raise ValueError("gene IDs not matching")

            # Add PROTEIN ID to the dictionary:
            mane_protein_ID = line[0]
            if ">" not in mane_protein_ID or "." not in mane_protein_ID :
                raise ValueError("> or . not in mane_protein_ID")
            else:
                mane_protein_ID = mane_protein_ID[1:]
                add_dict_branch(dict_MANE,[mane_gene_ID,"mane_protein_ID"],mane_protein_ID)

        # Add lines below the header (= cDNA) to the dictionary, skip for clinical IDs:
        else:
            if skip_line:
                continue
            dict_MANE[mane_gene_ID]["protein"] += line.rstrip()


###################################################################################################
# CREATE MANE GTF DICTIONARY ######################################################################
###################################################################################################

# Read MANE GTF file to get UTRs and CDS:
dict_MANE_positions = {}
dict_CDS_UTR_overlap = {}
dict_exon_overlap = {}

print("\n>>> READING MANE GTF FILE")
with open(data_path + "/MANE/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf", 'r') as gtffile:
    for counter, line in enumerate(gtffile):
        line = line.split()
        
        feature = line[2]
        # Skip lines with genes:
        if feature == "gene":
            continue
        
        # Define GTF file structure:
        gtf_start = int(line[3])
        gtf_end = int(line[4])
        gtf_score = line[5]
        gtf_strand = line[6]
        gtf_frame = line[7]
        mane_gene_ID = line[9].split('"')[1][:15]
        mane_transcript_ID = line[11].split('"')[1]
        mane_gene_name = line[15].split('"')[1]
        gtf_tags = line[16:]

        # Skip lines with clinical IDs:
        if '"MANE_Plus_Clinical";' in gtf_tags:
            continue
        
        # For transcipts, add their span to the dictionary:
        if feature == "transcript":
            add_dict_branch(dict_MANE_positions,[mane_gene_ID,"transcript_span"],range(gtf_start,gtf_end+1))
            add_dict_branch(dict_MANE_positions,[mane_gene_ID,"gtf_strand"],gtf_strand)
            continue

        exon_nr = int(line[21][:-1])
        add_dict_branch(dict_MANE_positions,[mane_gene_ID,"mane_gene_ID"],mane_gene_ID)
        add_dict_branch(dict_MANE_positions,[mane_gene_ID,"mane_transcript_ID"],mane_transcript_ID)

        if feature == "exon":
            add_dict_branch(dict_MANE_positions,[mane_gene_ID, feature],[[]])
            while len(dict_MANE_positions[mane_gene_ID][feature])<exon_nr:
                dict_MANE_positions[mane_gene_ID][feature].append([])
            dict_MANE_positions[mane_gene_ID][feature][exon_nr-1]=range(gtf_start,gtf_end+1)
        else:
            add_dict_branch(dict_MANE_positions,[mane_gene_ID, feature],[])
            dict_MANE_positions[mane_gene_ID][feature].append(range(gtf_start,gtf_end+1))


###################################################################################################
#  CDS AND UTR LENGTHS ############################################################################
###################################################################################################

for mane_gene_ID in dict_MANE_positions.keys():
    # Calculate UTR LENGTHS and add to dictionary:
    utr3 = 0
    utr5 = 0
    for pos in dict_MANE_positions[mane_gene_ID]["UTR"]:
        if dict_MANE_positions[mane_gene_ID]["gtf_strand"] == "+":
            if pos[0] < dict_MANE_positions[mane_gene_ID]["CDS"][0][0]:
                utr5 += len(pos)
                add_dict_branch(dict_MANE_positions,[mane_gene_ID,"UTR5"],[])
                dict_MANE_positions[mane_gene_ID]["UTR5"].append(pos)
            if pos[0] > dict_MANE_positions[mane_gene_ID]["CDS"][0][0]:
                utr3 += len(pos)
                add_dict_branch(dict_MANE_positions,[mane_gene_ID,"UTR3"],[])
                dict_MANE_positions[mane_gene_ID]["UTR3"].append(pos)
        if dict_MANE_positions[mane_gene_ID]["gtf_strand"] == "-":
            if pos[0] > dict_MANE_positions[mane_gene_ID]["CDS"][0][0]:
                utr5 += len(pos)
                add_dict_branch(dict_MANE_positions,[mane_gene_ID,"UTR5"],[])
                dict_MANE_positions[mane_gene_ID]["UTR5"].append(pos)
            if pos[0] < dict_MANE_positions[mane_gene_ID]["CDS"][0][0]:
                utr3 += len(pos)
                add_dict_branch(dict_MANE_positions,[mane_gene_ID,"UTR3"],[])
                dict_MANE_positions[mane_gene_ID]["UTR3"].append(pos)
    add_dict_branch(dict_MANE_positions,[mane_gene_ID,"UTR3_length"],utr3)
    add_dict_branch(dict_MANE_positions,[mane_gene_ID,"UTR5_length"],utr5)

    # Calculate CDS LENGTHS and add to dictionary:
    cds = 0
    for pos in dict_MANE_positions[mane_gene_ID]["CDS"]:
        cds += len(pos)
    add_dict_branch(dict_MANE_positions,[mane_gene_ID,"CDS_length"],cds)

    # Calculate INTRON LENGTHS and add to dictionary:
    transcript_length = len(dict_MANE_positions[mane_gene_ID]["transcript_span"])
    intron_length = transcript_length - dict_MANE_positions[mane_gene_ID]["UTR3_length"] - dict_MANE_positions[mane_gene_ID]["UTR5_length"] - dict_MANE_positions[mane_gene_ID]["CDS_length"]
    add_dict_branch(dict_MANE_positions,[mane_gene_ID,"intron_length"],intron_length)

    # Calculate cDNA LENGTHS and add to dictionary:
    cDNA_length = len(dict_MANE[mane_gene_ID]["cDNA"])
    add_dict_branch(dict_MANE_positions,[mane_gene_ID,"cDNA_length"],cDNA_length)


for gene in dict_MANE.keys():
    add_dict_branch(dict_MANE,[gene,"UTR5"],"")
    dict_MANE[gene]["UTR5"] = dict_MANE[gene]["cDNA"][:dict_MANE_positions[gene]["UTR5_length"]]

    add_dict_branch(dict_MANE,[gene,"UTR3"],"")
    dict_MANE[gene]["UTR3"] = dict_MANE[gene]["cDNA"][-dict_MANE_positions[gene]["UTR3_length"]:]

    add_dict_branch(dict_MANE,[gene,"CDS"],"")
    dict_MANE[gene]["CDS"] = dict_MANE[gene]["cDNA"][dict_MANE_positions[gene]["UTR5_length"]:-dict_MANE_positions[gene]["UTR3_length"]]


# %%
###################################################################################################
# TESTS ###########################################################################################
###################################################################################################
"""
print(mane_clinical_IDs)

print("\ncDNA is EMPTY in: ##################################")
for rna in dict_MANE.keys():
    if dict_MANE[rna]["cDNA"] == "":
        print(dict_MANE[rna]["mane_gene_ID"], dict_MANE[rna]["mane_gene_name"])

print("\nCDS is EMPTY in: ##################################")
for rna in dict_MANE.keys():
    if dict_MANE[rna]["CDS"] == "":
        print(dict_MANE[rna]["mane_gene_ID"], dict_MANE[rna]["mane_gene_name"])

my_length = 4
print(f"\nUTR3 is < {my_length} in: ##################################")
for rna in dict_MANE.keys():
    if len(dict_MANE[rna]["UTR3"]) < my_length:
        print(dict_MANE[rna]["mane_gene_ID"], dict_MANE[rna]["mane_gene_name"])

print("\nUTR5 is EMPTY in: ##################################")
for rna in dict_MANE.keys():
    if dict_MANE[rna]["UTR5"] == "":
        print(dict_MANE[rna]["mane_gene_ID"], dict_MANE[rna]["mane_gene_name"])

# %%
"""