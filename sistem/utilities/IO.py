import numpy as np
from sklearn.preprocessing import MinMaxScaler
from pyfaidx import Fasta
import pickle as pkl
import os
import logging

from sistem.utilities.utilities import init_mbit

def load_gs(gs_path):
    if not os.path.exists(gs_path):
        raise ValueError("No file exists at the provided path")
    with open(gs_path, 'rb') as f:
        gs = pkl.load(f)
    max_region_id = max(gs.anatomy.libraries[0].regions.values()) 
    init_mbit(max_region_id.bit_length())
    return gs

def load_distances_from_file(nsites, path):
    with open(path) as f:
        lines = [l.strip().split() for l in f.readlines()]
        if len(lines) == 1:
            dists = np.array(lines[0], dtype=float).reshape(-1, 1)
            norm_dists = MinMaxScaler().fit_transform(dists)
            return norm_dists
        else:
            idxs = np.triu_indices(nsites)
            dists = np.array(lines, dtype=float)[idxs]
            norm_dists = MinMaxScaler().fit_transform(dists)
            return norm_dists

def read_fasta(input_fasta):
    ref = Fasta(input_fasta, one_based_attributes=False)
    return ref

def setup_logger(file_path=None):
    if file_path:
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        filehandler = logging.FileHandler(file_path, 'w')
        filehandler.setFormatter(logging.Formatter('%(message)s'))
        for hdlr in logger.handlers[:]:
            if isinstance(hdlr,logging.FileHandler):
                logger.removeHandler(hdlr)
        logger.addHandler(filehandler)
    else:
        logging.basicConfig(level=logging.INFO, format='%(message)s')