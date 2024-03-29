import sys
import numpy as np
import pandas as pd
import os
import glob
from srr_functions import create_srr_to_cell

USAGE = "python make_rsem_dataframe.py <level> <gtf> <srrLookup> <seriesMatrix> <outCountsFile> <outTpmsFile> <namesOut> <idsOut>"
if len(sys.argv) != 9: print(USAGE); exit();
level, gtf, srr_lookup, series_matrix, outcounts, outtpms, outn, outids = sys.argv[1:]

print(f"getting gene names from {gtf}")
with open(gtf) as gtfhandle:
    line_ct = sum([1 for line in gtfhandle])
gene_id_to_name, gene_id_to_biotype = {}, {}
transcript_id_to_name, transcript_id_to_biotype, transcript_id_to_gene = {}, {}, {}
with open(gtf) as gtfhandle:
    for i, line in enumerate(gtfhandle):
        if i % 50000 == 0: print(f"Processed {i} lines out of {line_ct} from {gtf}.")
        if not line.startswith('#'):
            line = line.strip().rstrip(';').split('\t')
            if len(line) < 9 or line[2] != "transcript": continue
            attribute_list = line[8].replace('; ', ';').split(';')
            attributes = {}
            for item in attribute_list:
                item = [x.strip('"') for x in item.replace(" \"",";").split(';')]
                if len(item) < 2:
                    print(f"Error: not a key-value pair: {attributes} {item} in line:\n{line}")
                attributes[item[0]] = item[1]
            if not "gene_name" in attributes:
                attributes["gene_name"] = attributes["gene_id"]
            if not "transcript_name" in attributes:
                attributes["transcript_name"] = attributes["transcript_id"]
            if any([x not in attributes for x in ["gene_id", "gene_name", "transcript_id", "transcript_name"]]):
                print("Error: \"gene_id\", \"gene_name\", \"transcript_id\", or \"transcript_name\" not found in line:\n" + "\t".join(line))
                exit(1)
            gene_id_to_name[attributes["gene_id"].strip("gene:")] = attributes["gene_name"]
            gene_id_to_biotype[attributes["gene_id"].strip("gene:")] = attributes["gene_biotype"] if "gene_biotype" in attributes else ""
            transcript_id_to_name[attributes["transcript_id"].strip("transcript:")] = attributes["transcript_name"]
            transcript_id_to_biotype[attributes["transcript_id"].strip("transcript:")] = attributes["transcript_biotype"] if "transcript_biotype" in attributes else ""
            transcript_id_to_gene[attributes["transcript_id"].strip("transcript:")] = attributes["gene_id"]

print(f"globbing ../results/quant/*.{level}.results")
files = glob.glob(f"../results/quant/*.{level}.results")
line_ct = sum(1 for line in open(files[0]))

print("Making SRR-to-cell lookup...")
srr_to_cell = create_srr_to_cell(series_matrix, srr_lookup)

def get_prefix(file):
    srr = os.path.basename(file).split(".")[0].split("_")[0]
    return srr_to_cell[srr]

doUseGene = level.startswith("gene")
ids = [line.split('\t')[0].strip("gene:").strip("transcript:") for line in open(files[0])]
ids[0] = "gene_id" if doUseGene else "transcript_id"
names = [
    (gene_id_to_name[id] if id in gene_id_to_name else id) if doUseGene else \
    (transcript_id_to_name[id] if id in transcript_id_to_name else id) \
    for id in ids[1:]]
biotypes = [
    (gene_id_to_biotype[id] if id in gene_id_to_biotype else "") if doUseGene else \
    (transcript_id_to_biotype[id] if id in transcript_id_to_biotype else "") \
    for id in ids[1:]]

def save_output(outfilename, id_array, in_files, col_idx, doUseGene):
    print(f"Reading data files to create {outfilename} ...")

    with open(outfilename, "w") as ensg_output_file, open(outfilename + ".ercc.csv", "w") as ercc_output_file:
        ensg_output_file.write(",".join([""] + id_array[1:]) + "\n")
        ercc_output_file.write(",".join([""] + id_array[1:]) + "\n")
        prefixes = [get_prefix(file) for file in in_files]
        file_order = np.argsort(prefixes)
        for file in np.array(in_files)[file_order]:
            prefix = get_prefix(file)
            print(f"reading {prefix}")

            df = pd.read_csv(file, sep='\t', header=None, usecols=[0, col_idx], names=['id', 'value'])
            df['id'] = df['id'].str.strip('gene:').str.strip('transcript:')
            df = df.set_index('id').transpose()
            df.index = [prefix]

            ensg_df = df.filter(regex="ENSG" if doUseGene else "ENST")
            ercc_df = df.filter(regex="ERCC")

            ensg_output_file.write(ensg_df.to_csv(header=False, index=True))
            ercc_output_file.write(ercc_df.to_csv(header=False, index=True))

    print(f"Data saved to {outfilename} and {outfilename}.ercc.csv")

for (filename, col_idx) in [(outcounts, 4), (outtpms, 5)]:
    save_output(filename, ids, files, col_idx, doUseGene)

print(f"Saving to {outn} ...")
np.savetxt(outn, np.column_stack((ids[1:], names, biotypes)), delimiter=",", fmt="%s")

print(f"Saving to {outids} ...")
np.savetxt(outids, 
    np.column_stack((list(transcript_id_to_gene.keys()), 
        list(transcript_id_to_gene.values()))), delimiter=",", fmt="%s")

print("Done.")
