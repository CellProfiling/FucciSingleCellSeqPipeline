import scvelo as scv
import sys, io
import numpy as np
from srr_functions import create_srr_to_cell

aloomfile = sys.argv[1]
srr_lookup = sys.argv[2]
series_matrix = sys.argv[3]
output = sys.argv[4]

aloom = scv.read_loom(aloomfile)
aobs = aloom.obs_names

srr_to_cell = create_srr_to_cell(series_matrix, srr_lookup)

# Use aobs to look up cell name from cells using the SRR link somehow
acells = [srr_to_cell[x.split(':')[1].split('Aligned')[0]] for x in aobs]
with open(output, 'w') as output:
    output.write('well_plate\n')
    output.write('\n'.join(acells))
