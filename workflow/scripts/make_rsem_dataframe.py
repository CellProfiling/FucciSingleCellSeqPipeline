import sys
import numpy as np
import pandas as pd
import os
import glob

USAGE = "python make_rsem_dataframe.py <level> <gff> <srrLookup> <seriesMatrix> <outCountsFile> <outTpmsFile> <namesOut> <idsOut>"
if len(sys.argv) != 9: print(USAGE); exit();
level, gff, srr_lookup, series_matrix, outcounts, outtpms, outn, outids = sys.argv[1:]

print(f"getting gene names from {gff}")
with open(gff) as gffhandle:
    line_ct = sum([1 for line in gffhandle])
gene_id_to_name, gene_id_to_biotype, gene_id_to_description = {}, {}, {}
transcript_id_to_name, transcript_id_to_biotype, transcript_id_to_description, transcript_id_to_gene = {}, {}, {}, {}
with open(gff) as gffhandle:
    for i, line in enumerate(gffhandle):
        if i % 50000 == 0: print(f"Processed {i} lines out of {line_ct} from {gff}.")
        if not line.startswith('#'):
            line = line.strip().rstrip(';').split('\t')
            if len(line) < 9: continue
            attribute_list = line[8].replace('; ', ';').split(';')
            attributes = {}
            for item in attribute_list:
                item = item.split('=')
                attributes[item[0]] = item[1]
            if "ID" in attributes and "Name" in attributes and attributes["ID"].startswith("gene:"):
                gene_id_to_name[attributes["ID"][len("gene:"):]] = attributes["Name"]
                gene_id_to_biotype[attributes["ID"][len("gene:"):]] = attributes["biotype"] if "ID" in attributes and "biotype" in attributes else ""
                gene_id_to_description[attributes["ID"][len("gene:"):]] = attributes["description"] if "ID" in attributes and "description" in attributes else ""
            elif "ID" in attributes and "Name" in attributes and attributes["ID"].startswith("transcript:"):
                transcript_id_to_name[attributes["ID"][len("transcript:"):]] = attributes["Name"]
                transcript_id_to_biotype[attributes["ID"][len("transcript:"):]] = attributes["biotype"] if "ID" in attributes and "biotype" in attributes else ""
                transcript_id_to_description[attributes["ID"][len("transcript:"):]] = attributes["description"] if "ID" in attributes and "description" in attributes else ""
                transcript_id_to_gene[attributes["ID"][len("transcript:"):]] = attributes["Parent"][len("gene:"):] if "ID" in attributes and "Parent" in attributes else ""



print(f"globbing ../results/quant/*.{level}.results")
files = glob.glob(f"../results/quant/*.{level}.results")
line_ct = sum(1 for line in open(files[0]))

print("Making SRR-to-cell lookup...")
srr_lookup = dict([x.split()[::-1] for x in open(srr_lookup)])
series_matrix_lines = open(series_matrix).readlines()
cells, srr, has_srr = [], [], []
for line in series_matrix_lines:
    if line.startswith('!Sample_title'):
        cells=[x.strip('"').strip('Single U2OS cell ').strip() for x in line.split('\t')[1:]]
    if line.startswith('!Sample_relation') and "SRX" in line:
        srrs = [x.split('=')[1][:-1].strip() for x in line.split('\t')[1:]]
        has_srr = [x in srr_lookup for x in srrs]
        srr = [srr_lookup[x] for x in srrs if x in srr_lookup]
cells=np.array(cells)[has_srr] # for testing
if len(cells) != len(srr):
    print("Error: series matrix has different cell and srr lengths")
    exit(1)
srr_to_cell = dict((srr[idx], cells[idx]) for idx in range(len(srr)))

def get_prefix(file):
    srr = os.path.basename(file).split(".")[0].split("_")[0]
    return srr_to_cell[srr]

doUseGene = level.startswith("gene")
ids = [line.split('\t')[0].strip("gene:").strip("transcript:") for line in open(files[0])]
ids[0] = "gene_id" if doUseGene else "transcript_id"
names = [gene_id_to_name[id] if doUseGene else transcript_id_to_name[id] for id in ids[1:]]
biotypes = [gene_id_to_biotype[id] if doUseGene else transcript_id_to_biotype[id] for id in ids[1:]]
description = [gene_id_to_description[id] if doUseGene else transcript_id_to_description[id] for id in ids[1:]]
counts = [np.array(ids)]
tpms = [np.array(ids)]
for file in files:
    counts.append(np.array([line.split('\t')[4] for line in open(file)]))
    tpms.append(np.array([line.split('\t')[5] for line in open(file)]))
    print("reading "+ get_prefix(file))

print(f"Saving to {outcounts} ...")
dataframe = np.row_stack(counts)
dataframe[1:,0] = [get_prefix(file) for file in files]
pddf = pd.DataFrame(dataframe[1:,1:], index=dataframe[1:,0], columns=dataframe[0,1:]).sort_index()
pddf.filter(regex="ENSG" if doUseGene else "ENST").to_csv(outcounts)
pddf.filter(regex="ERCC").to_csv(outcounts + ".ercc.csv")

print(f"Saving to {outtpms} ...")
dataframe = np.row_stack(tpms)
dataframe[1:,0] = [get_prefix(file) for file in files]
pddf = pd.DataFrame(dataframe[1:,1:], index=dataframe[1:,0], columns=dataframe[0,1:]).sort_index()
pddf.filter(regex="ENSG" if doUseGene else "ENST").to_csv(outtpms)
pddf.filter(regex="ERCC").to_csv(outtpms + ".ercc.csv")

print(f"Saving to {outn} ...")
np.savetxt(outn, np.column_stack((ids[1:], names, biotypes, description)), delimiter=",", fmt="%s")

print(f"Saving to {outids} ...")
np.savetxt(outids, np.column_stack((list(transcript_id_to_gene.keys()), list(transcript_id_to_gene.values()))), delimiter=",", fmt="%s")

print("Done.")
