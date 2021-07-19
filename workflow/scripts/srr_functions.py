import re
import numpy as np

def create_srr_to_cell(series_matrix_file, srr_lookup_file):
    srr_lookup = dict([x.split()[::-1] for x in open(srr_lookup_file)])
    series_matrix_lines = open(series_matrix_file).readlines()
    cells, srr, has_srr = [], [], []
    for line in series_matrix_lines:
        if line.startswith('!Sample_title'):
            pattern = re.compile(r'cell ([A-Z][0-9]+_[0-9]+)')
            cells = [pattern.search(x).group(1) for x in line.split('\t')[1:]]
        if line.startswith('!Sample_relation') and "SRX" in line:
            pattern = re.compile(r'=(SRX[0-9]+)')
            srrs = [pattern.search(x).group(1) for x in line.split('\t')[1:]]
            has_srr = [x in srr_lookup for x in srrs]
            srr = [srr_lookup[x] for x in srrs if x in srr_lookup]
    cells=np.array(cells)[has_srr] # for testing
    if len(cells) != len(srr):
        raise RuntimeError("Error: series matrix has different cell and srr lengths")

    srr_to_cell = dict((srr[idx], cells[idx]) for idx in range(len(srr)))
    return srr_to_cell
