import numpy as np
import os

from sistem.data.utils import bin_region_values, get_mutated_basepairs, count_region_CNs
from sistem.lineage import CNA, SNV
from sistem.utilities.utilities import get_reg_id, get_bin_coordinates, sort_chrom_names
from sistem.utilities.IO import read_fasta

def save_mutations_from_tree(tree, out_dir, bin_size=None):
    log_mutations(tree, out_dir)
    save_CNPs(tree, out_dir, bin_size=bin_size)
    save_SNVs(tree, out_dir)

def save_CNPs(tree, out_dir, bin_size):
    """
    Saves the CNPs for the cells in a tree. Stores in two separate files for observed and ancestral cells. If region_len != bin_size, saves the region CNPs in a separate file.
    """
    chrom_lens = tree.root.library.chrom_lens
    region_len = tree.root.library.region_len
    chrnames = sort_chrom_names(chrom_lens.keys())

    if bin_size == None:
        bin_size = region_len

    region_coords = get_bin_coordinates(chrom_lens, region_len)

    observed_cells = [node for node in tree.iter_leaves()]
    ancestral_cells = [node for node in tree.iter_ancestors()]
    
    for group, keyword in zip([observed_cells, ancestral_cells], ['observed_', 'ancestral_']):
        region_profiles = {}
        for cell in group:
            region_profiles[cell.name] = count_region_CNs(cell)
        
        if bin_size == region_len:
            CNPs_to_tsv(region_profiles, region_coords, chrnames, os.path.join(out_dir, f'{keyword}CNPs.tsv'))
        else:
            #CNPs_to_tsv(region_profiles, region_coords, chrnames, os.path.join(out_dir, f'region_{keyword}CNPs.tsv'))

            bin_profiles = {}
            for cell in group:
                bin_profiles[cell.name] = bin_region_CNs(region_profiles[cell.name], region_len, bin_size)
            bin_coords = get_bin_coordinates(chrom_lens, bin_size)
            CNPs_to_tsv(bin_profiles, bin_coords, chrnames, os.path.join(out_dir, f'{keyword}CNPs.tsv'))

def CNPs_to_tsv(profiles, coords, chrnames, out_path):
    cell_names = list(profiles.keys())
    cell_names.sort()

    headers = ['Cell/Clone', 'Chrom', 'Start', 'End', 'Hap CN', 'Total CN']
    rows = []
    for chrname in chrnames:
        for i, (start, end) in enumerate(coords[chrname]):
            for cell in cell_names:
                CNa, CNb = round(profiles[cell][chrname][0][i], 3), round(profiles[cell][chrname][1][i], 3)
                CNtot = round(CNa + CNb, 3)
                rows.append(f'{cell}\t{chrname}\t{start}\t{end}\t{CNa},{CNb}\t{CNtot}\n')
    
    with open(out_path, 'w+') as f:
        f.write('\t'.join(headers) + '\n')
        for row in rows:
            f.write(row)

def bin_region_CNs(CNs, region_len, bin_size, total=False):
    """
    Partitions region CNs into continuous bins with the resulting bin CN being the average copy number of all regions belonging to that bin. The ratio bin_size/region_len may be a float value, in which case a weighted average is used and continuity of regions is maintained.

    Args:
        CNs (dict): Dictionary containing region copy numbers.
        region_len (int): Length of regions.
        bin_size (int): Length of bins.
    
    Returns:
        bin_CNs (dict): Binned copy numbers.
    """
    bin_CNs = {}
    r2b = bin_size/region_len

    for chrname in CNs:
        if total:
            bin_CNs[chrname] = bin_region_values(CNs[chrname], r2b)
        else:
            bin_CNs[chrname] = [bin_region_values(CNs[chrname][0], r2b), bin_region_values(CNs[chrname][1], r2b)]
    
    return bin_CNs


def save_SNVs(tree, out_dir, ref = None, alt_ref = None):
    """
    ref (str, None): Path to reference genome
    alt_ref (str, None):: Path to alt reference genome
    """
    chrom_lens = tree.root.library.chrom_lens
    region_len = tree.root.library.region_len
    chrnames = sort_chrom_names(chrom_lens.keys())

    observed_cells = [node for node in tree.iter_leaves()]
    mutated_bps = get_mutated_basepairs(observed_cells, chrnames, region_len)

    ref_handler, alt_ref_handler = None, None
    if ref:
        ref_handler = read_fasta(ref)
    if alt_ref:
        alt_ref_handler = read_fasta(alt_ref)

    profiles = create_SNV_profiles(observed_cells, mutated_bps, region_len, ref_handler=ref_handler, alt_ref_handler=alt_ref_handler)

    with open(os.path.join(out_dir, 'SNV_profiles.tsv'), 'w+') as f:
        f.write(f"Cell/Clone\tChrom\tPos\tGT\n")
        for i,(chrname, pos) in enumerate(mutated_bps):
            for cell in observed_cells:
                x = profiles[cell][i]
                f.write(f"{cell.name}\t{chrname}\t{pos}\t{x[0]}/{x[1]}\n")

def create_SNV_profiles(cells, mutated_bps, region_len, ref_handler = None, alt_ref_handler = None):
    n = len(mutated_bps)
    nucs = ['A', 'C', 'G', 'T']
    if ref_handler is None and alt_ref_handler is None:
        A_nucs = [str(np.random.choice(nucs)) for _ in range(n)]
        B_nucs = A_nucs
    else:
        A_nucs = []
        for chrname, pos in mutated_bps:
            bp = ref_handler[chrname][pos].seq.upper()
            if bp == 'N':
                A_nucs.append(str(np.random.choice(nucs)))
            else:
                A_nucs.append(bp)
        if alt_ref_handler is None:
            B_nucs = A_nucs
        else:
            B_nucs = []
            for i,(chrname, pos) in enumerate(mutated_bps):
                bp = alt_ref_handler[chrname][pos].seq.upper()
                if bp == 'N':
                    B_nucs.append(A_nucs[i])
                else:
                    B_nucs.append(bp)
    ref_nucs = [A_nucs, B_nucs]

    profiles = {cell: [['.', '.'] for _ in range(n)] for cell in cells}
    for cell in cells:
        for i,(chrname, pos) in enumerate(mutated_bps):
            r = pos // region_len
            p = pos % region_len
            for chromosome in cell.genome[chrname]:
                if profiles[cell][i][chromosome.allele] == '.':
                    for q in chromosome.seq:
                        if get_reg_id(q) == r:
                            profiles[cell][i][chromosome.allele] = ref_nucs[chromosome.allele][i]
                            break

                for q,positions in chromosome.SNVs.items():
                    if get_reg_id(q) == r:
                        if p in positions:
                            alt_nucs = nucs[:]
                            alt_nucs.remove(ref_nucs[chromosome.allele][i])
                            alt = alt_nucs[positions[p]]
                            profiles[cell][i][chromosome.allele] = alt
    return profiles

def log_mutations(tree, out_dir):
    with open(os.path.join(out_dir, 'CNA_events.tsv'), 'w+') as f:
        f.write(f'Cell/Clone\tType\tDriver\tChrom\tAllele\tArm\tStart Index\tLength\tCopies\tRef Region Indices\n')
        for node in tree.iter_descendants():
            for e in node.events:
                if isinstance(e, CNA):
                    if e.category == 'focal_amplification' or e.category == 'focal_deletion':
                        f.write(f"{node.name}\t{e.category}\t{e.driver}\t{e.chrom}\t{e.allele}\t{e.arm}\t{e.start}\t{e.end - e.start}\t{e.copies}\t{','.join([str(get_reg_id(g)) for g in e.gene_ids])}\n")
                    else:
                        f.write(f"{node.name}\t{e.category}\t{e.driver}\t{e.chrom}\t{e.allele}\t{e.arm}\tNone\tNone\tNone\tNone\n")
    
    with open(os.path.join(out_dir, 'SNV_events.tsv'), 'w+') as f:
        f.write(f'Cell/Clone\tDriver\tChrom\tAllele\tRegion Index\tRef Region Index\tPosition\tBP\n')
        for node in tree.iter_descendants():
            for e in node.events:
                if isinstance(e, SNV):
                    f.write(f"{node.name}\t{e.driver}\t{e.chrom}\t{e.allele}\t{e.index}\t{get_reg_id(e.region)}\t{e.position}\t{e.bp}\n")

def save_singlecell_readcounts(readcounts, out_path, chrom_lens, bin_size):
    chrnames = sort_chrom_names(chrom_lens.keys())
    bin_coords = get_bin_coordinates(chrom_lens, bin_size)
    with open(os.path.join(out_path, 'readcounts.tsv'), 'w+') as f:
        headers = ['Cell', 'Chrom', 'Start', 'End', 'Acount', 'Bcount']
        f.write('\t'.join(headers) + '\n')

        for chrname in chrnames:
            for i, (start, end) in enumerate(bin_coords[chrname]):
                for cell, rcs in readcounts.items():
                    f.write(f"{cell.name}\t{chrname}\t{start}\t{end}\t{rcs[chrname][0][i]}\t{rcs[chrname][1][i]}\n")

def save_clonal_readcounts(mutated_basepairs, readcounts, site_ids, out_path):
    with open(os.path.join(out_path, 'readcounts.tsv'), 'w+') as f:
        headers = ['Chrom', 'Pos'] + site_ids
        f.write('\t'.join(headers) + '\n')

        for i,(chrname,pos) in enumerate(mutated_basepairs):
            counts = [f"{readcounts[s][i][0]},{readcounts[s][i][1]}" for s in range(len(site_ids))]
            counts_combined = '\t'.join(counts)
            f.write(f"{chrname}\t{pos}\t{counts_combined}\n")

def site_CN_averages(tree, out_dir, bin_size):
    chrom_lens = tree.root.library.chrom_lens
    region_len = tree.root.library.region_len
    chrnames = sort_chrom_names(chrom_lens.keys())

    bin_coords = get_bin_coordinates(chrom_lens, bin_size)

    observed_clones = [node for node in tree.iter_leaves()]
    nsites = len(set([clone.site for clone in observed_clones]))
    
    for clone in observed_clones:
        if 'count' not in clone.info:
            raise ValueError("Leaves in Tree have not been initialized with a 'count' key/pair in the info attribute. Ensure you have generated the tree with GrowthSimulator.simulate_clonal_lineage, as this function does not work on trees generated with GrowthSimulator.simulate_singlecell_lineage.")

    cell_counts = {}
    for s in range(nsites):
        cell_counts[s] = {}
        for clone in observed_clones:
            if clone.site == s:
                cell_counts[s][clone] = clone.info['count']

    #Add diploid fractions to metastatic sites
    dip = tree.find('diploid')
    if dip is not None:
        dip_count = cell_counts[0][dip]
        prim_tumor_cells = sum(list(cell_counts[0].values())) - dip_count
        ratio = dip_count / prim_tumor_cells
        for s in range(1, nsites):
            met_tumor_cells = sum(list(cell_counts[s].values()))
            if met_tumor_cells > 0:
                num_dip = max(round(ratio*met_tumor_cells), 1)
                cell_counts[s][dip] = num_dip
    
    CN_averages = {s: {} for s in range(nsites)}
    for s in range(nsites):
        cur_clones = list(cell_counts[s].keys())
        weights = [cell_counts[s][clone] for clone in cur_clones]

        profiles = {}
        for clone in cur_clones:
            region_profiles = count_region_CNs(clone, total=True)
            profiles[clone] = bin_region_CNs(region_profiles, region_len, bin_size, total=True)
        
        for chrname in chrnames:
            CNs_flat = np.array([profiles[clone][chrname] for clone in cur_clones])
            chr_CN_avgs = np.average(CNs_flat, axis=0, weights=weights)
            CN_averages[s][chrname] = chr_CN_avgs
    
    site_ids = []
    for s in range(nsites):
        cur_clones = [clone.name for clone in cell_counts[s].keys()]
        if 'diploid' in cur_clones:
            cur_clones.remove('diploid')
        name = cur_clones[0]
        site_id = ''
        for char in name[:2]:
            if char.isupper():
                site_id += char
        site_ids.append(site_id)

    with open(os.path.join(out_dir, 'CN_averages.tsv'), 'w+') as f:
        f.write(f"Chrom\tStart\tEnd\tSite\tCN\n")
        for chrname in chrnames:
            for i, (start, end) in enumerate(bin_coords[chrname]):
                for s in range(nsites):
                    f.write(f"{chrname}\t{start}\t{end}\t{s}\t{CN_averages[s][chrname][i]}\n")