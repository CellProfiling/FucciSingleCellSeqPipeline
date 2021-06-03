import scvelo as scv
import sys, io
import numpy as np

aloomfile = sys.argv[1]
srr_lookup = sys.argv[2]
series_matrix = sys.argv[3]
output = sys.argv[4]

aloom = scv.read_loom(aloomfile)
aobs = aloom.obs_names

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

# Use aobs to look up cell name from cells using the SRR link somehow
srr_to_cell = dict((srr[idx], cells[idx]) for idx in range(len(srr)))
acells = [srr_to_cell[x.split(':')[1].split('Aligned')[0]] for x in aobs]
with open(output, 'w') as output:
    output.write('well_plate\n')
    output.write('\n'.join(acells))
