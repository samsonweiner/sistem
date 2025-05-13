import numpy as np
import sys
import itertools
import shutil
from typing import Optional, List

class DependencyError(Exception):
    """Exception raised when required dependencies are missing."""
    pass

def check_dependencies(commands: Optional[List[str]] = None):
    missing = []
    for cmd in commands:
        if shutil.which(cmd) is None:
            missing.append(cmd)
    if missing:
        missing_combined = ', '.join(missing)
        raise DependencyError(f"The following binaries must be installed and callable: {missing_combined}.")
    return True

def compute_growth_rate(N0, Nt, tgen):
    '''
    Computes the exponential growth rate.

    Args:
        N0 (int): Population size at time t0.
        Nt (int): Population size at time t.
        tgen (int): Number of generations between t0 and t.
    
    Returns:
        (float): Exponential growth rate.
    '''
    return (np.log(Nt) - np.log(N0)) / tgen
    
def sort_chrom_names(chrnames):
    all_contains_prefix = np.all(np.array([chrname[:3] for chrname in chrnames]) == 'chr')
    all_missing_prefix = np.all(np.array([chrname[:3] for chrname in chrnames]) != 'chr')
    # If some have prefix 'chr' and others do not, naming is inconsistent and return what was inputted.
    if not (all_contains_prefix or all_missing_prefix):
        return chrnames
    
    if all_contains_prefix:
        # If any chrnames are are just 'chr', then return.
        if not np.all(np.array([len(chrname) for chrname in chrnames]) > 3):
            return chrnames
        chrnames = [chrname[3:] for chrname in chrnames]
    
    int_chrnames = [int(x) for x in chrnames if x.isdigit()]
    XY_chrnames = [x for x in chrnames if x[0] == 'X' or x[0] == 'Y']
    alt_chrnames = [x for x in chrnames if not (x.isdigit() or x[0] == 'X' or x[0] == 'Y')]
    int_chrnames.sort()
    XY_chrnames.sort()
    alt_chrnames.sort()
    sorted_chrnames = int_chrnames + XY_chrnames + alt_chrnames
    if all_contains_prefix:
        sorted_chrnames = [f'chr{x}' for x in sorted_chrnames]
    return sorted_chrnames
    
def get_bin_coordinates(chrom_lens, bin_size):
    bin_coords = {chrname: [] for chrname in chrom_lens}
    for chrname, chrlen in chrom_lens.items():
        nbins = int(chrlen // bin_size)
        rem = chrlen % bin_size
        if rem > 0:
            nbins += 1
        for i in range(nbins):
            start = int(i*bin_size + 1)
            end = int(min((i+1)*bin_size, chrlen))
            bin_coords[chrname].append((start, end))
    return bin_coords
    
def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


def combine_dicts(*dicts):
    """
    Given a set of 2-depth dictionaries of sets, combines all keys and subkeys into a single dictionary with unique values present in any of the inputs.

    Args:
        *dicts (iterable): Collection of dictionaries.
    
    Returns:
        combined (dict): Output dictionary.
    """

    combined = {}
    for d in dicts:
        for key in d:
            if key not in combined:
                combined[key] = {}
            for subkey in d[key]:
                if subkey not in combined[key]:
                    combined[key][subkey] = set()
                combined[key][subkey].update(d[key][subkey])

    #for key, subkeys in combined.items():
    #    for subkey in subkeys:
    #        combined[key][subkey] = list(set(combined[key][subkey]))
    
    return combined

def iter_by_chunk(iterable, chunksize):
    return itertools.zip_longest(*[iter(iterable)] * chunksize)

# Initialize the MBIT to be one significant bit higher than what is needed to represent the largest region id. All bits to the left of MBIT will be used as a unique id.
def init_mbit(val):
    global MBIT
    val = val << 1
    MBIT = val

# Given any integer x, return the unique id.
def get_reg_count(x):
    mask = ~((1 << (MBIT + 1)) - 1)
    y = (x & mask) >> (MBIT + 1)
    return y

# Given any integer x, return the true region id
def get_reg_id(x):
    mask = ((1 << (MBIT + 1)) - 1)
    r = (x & mask)
    return r

# Given true region id x and a unique identifier b, return x with b prepended to it before the MBIT.
def update(x, b):
    mask = ((1 << (MBIT + 1)) - 1)
    x_prime = (b << (MBIT + 1)) | (x & mask)
    return x_prime

# Given integer x, return x with MBIT flipped
#def inc(x):
#    return update(x, 1)