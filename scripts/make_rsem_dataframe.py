import sys
import numpy as np
import os

USAGE = "python make_rsem_dataframe.py <column_index> <outFile> file1 file2 ... fileN"

if len(sys.argv) < 4:
    print(USAGE)
    exit()

column = int(sys.argv[1])
outf, files = sys.argv[2], sys.argv[3:]
line_ct = sum(1 for line in open(files[0]))

ids = np.array([line.split('\t')[0].strip("gene:") for line in open(files[0])])
columns = [ids]
for file in files:
    columns.append(np.array([line.split('\t')[column] for line in open(file)]))
dataframe = np.column_stack(columns)
dataframe[0,:] = ["id"].extend([os.path.basename(file) + "_" + file.strip("_R1.fastq.gz") for file in files])

np.savetxt(outf, dataframe, delimiter="\t", fmt="%s")
