import re
import numpy as np

def load_srr_lookup(srr_lookup_file):
    """
    Load the SRR lookup file into a dictionary.
    
    Args:
        srr_lookup_file (str): Path to the SRR lookup file.
    
    Returns:
        dict: A dictionary containing SRR IDs as keys and corresponding values.
    """
    with open(srr_lookup_file) as file:
        srr_lookup = dict(line.strip().split()[::-1] for line in file)
    return srr_lookup

def parse_series_matrix_file(series_matrix_file, srr_lookup):
    """
    Parse the series matrix file to extract cell information and SRR IDs.
    
    Args:
        series_matrix_file (str): Path to the series matrix file.
        srr_lookup (dict): A dictionary containing SRR IDs as keys and corresponding values.
    
    Returns:
        tuple: A tuple containing the list of cells, the list of SRR IDs, and a boolean list
               indicating whether each SRR ID is in the srr_lookup dictionary.
    """
    with open(series_matrix_file) as file:
        series_matrix_lines = file.readlines()

    cells, srrs, has_srr = [], [], []
    cell_pattern = re.compile(r'cell ([A-Z][0-9]+_[0-9]+)')
    srr_pattern = re.compile(r'=(SRX[0-9]+)')

    for line in series_matrix_lines:
        if line.startswith('!Sample_title'):
            cells = [cell_pattern.search(x).group(1) for x in line.split('\t')[1:]]
        if line.startswith('!Sample_relation') and "SRX" in line:
            srrs = [srr_pattern.search(x).group(1) for x in line.split('\t')[1:]]
            has_srr = [srr in srr_lookup for srr in srrs]
            srr = [srr_lookup[srr] for srr in srrs if srr in srr_lookup]

    return cells, srrs, has_srr

def create_srr_to_cell(series_matrix_file, srr_lookup_file):
    """
    Create a dictionary that maps SRR IDs to cell information.
    
    Args:
        series_matrix_file (str): Path to the series matrix file.
        srr_lookup_file (str): Path to the SRR lookup file.
    
    Returns:
        dict: A dictionary that maps SRR IDs to cell information.
    """
    srr_lookup = load_srr_lookup(srr_lookup_file)
    cells, srrs, has_srr = parse_series_matrix_file(series_matrix_file, srr_lookup)

    cells = np.array(cells)[has_srr]

    if len(cells) != len(srrs):
        raise RuntimeError("Error: series matrix has different cell and srr lengths")

    srr_to_cell = {srr: cell for srr, cell in zip(srrs, cells)}
    return srr_to_cell
