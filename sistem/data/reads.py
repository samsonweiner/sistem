import os
import subprocess
import numpy as np
import multiprocessing
from Bio import SeqIO
from scipy.stats import poisson
from typing import Optional

from sistem.utilities.utilities import iter_by_chunk, get_reg_id, sort_chrom_names, check_dependencies
from sistem.data.utils import gen_coverage, get_alpha_beta, bin_region_readcounts, get_mutated_basepairs, count_region_CNs
from sistem.data.sequence import build_cell_ref
from sistem.data.profiles import save_singlecell_readcounts, save_clonal_readcounts
from sistem.lineage import Tree
from sistem.parameters import Parameters, fill_params
from sistem.utilities.IO import read_fasta
from sistem.utilities.utilities import init_mbit

# Draw readcounts for each bin
def draw_readcounts_nonuniform(num_windows, Aa, Bb, mean_rc, interval=5):
    cov_scales = gen_coverage(num_windows, interval, Aa, Bb)
    readcounts = [round(poisson.rvs(mean_rc*x)) for x in cov_scales]
    return readcounts

def gen_reads_cell(
    cell, 
    readcounts, 
    out_dir, 
    ref_path, 
    alt_ref_path, 
    read_len, 
    seq_error
):
    ref = read_fasta(ref_path)
    if alt_ref_path is None:
        refs = [ref, ref]
    else:
        alt_ref = read_fasta(alt_ref_path)
        refs = [ref, alt_ref]

    chrom_lens = cell.library.chrom_lens
    regions = cell.library.regions
    region_len = cell.library.region_len

    max_region_id = max(regions.values())
    init_mbit(max_region_id.bit_length())

    last_region_lens = {chrname: int(l % region_len) for chrname,l in chrom_lens.items()}

    read_dir = os.path.join(out_dir, 'reads')
    prefix = os.path.join(read_dir, f"{cell.name}")
    prefixes = [f"{prefix}_allele0", f"{prefix}_allele1"]
    for cur_prefix in prefixes:
        open(f"{cur_prefix}.read1.fastq.gz", 'w+').close()
        open(f"{cur_prefix}.read2.fastq.gz", 'w+').close()

    for chrname in chrom_lens:
        for chromosome in cell.genome[chrname]:
            if len(chromosome) > 0:
                cur_ref = refs[chromosome.allele]
                cur_prefix = prefixes[chromosome.allele]
                cur_readcounts = readcounts[chromosome]

                chrom_fa_path = f"{cur_prefix}.{chrname}_{chromosome.homolog_id}.fa"
                build_cell_ref(chromosome, cur_ref, regions[chrname], region_len, chrom_fa_path)

                start,end = 0,0
                for i,r in enumerate(chromosome.seq):
                    q = get_reg_id(r)
                    rc = cur_readcounts[i]
                    if q == regions[chrname] - 1:
                        end += last_region_lens[chrname]
                    else:
                        end += region_len
                    
                    region_fa_path = f"{chrom_fa_path[:-3]}_{int(start)}_{int(end)}.fa"
                    with open(region_fa_path, 'w+') as f:
                        call = subprocess.run(['samtools', 'faidx', chrom_fa_path, f"{chrname}:{int(start)}-{int(end)}"], stdout=f)

                    for record in SeqIO.parse(region_fa_path, 'fasta'):
                        total_N = record.seq.count('N')
                    total_bp = end - start + 1
                    N_ratio = total_N / total_bp
                    rc = int(rc*(1-N_ratio))

                    if rc > 0:
                        call = subprocess.run(['dwgsim', '-H', '-o', '1', '-N', str(rc), '-1', str(read_len), '-2', str(read_len), '-e', str(seq_error), '-E', str(seq_error), region_fa_path, region_fa_path[:-3]], capture_output=True, text=True)

                        with open(f"{cur_prefix}.read1.fastq.gz", 'a') as f1, open(f"{cur_prefix}.read2.fastq.gz", 'a') as f2:
                            call = subprocess.run(['cat', region_fa_path[:-3] + '.bwa.read1.fastq.gz'], stdout=f1)
                            call = subprocess.run(['cat', region_fa_path[:-3] + '.bwa.read2.fastq.gz'], stdout=f2)
                        os.remove(region_fa_path[:-3] + '.bwa.read1.fastq.gz')
                        os.remove(region_fa_path[:-3] + '.bwa.read2.fastq.gz')
                        os.remove(region_fa_path[:-3] + '.mutations.txt')
                        os.remove(region_fa_path[:-3] + '.mutations.vcf')
                    os.remove(region_fa_path)
                    start = end
                os.remove(chrom_fa_path)
                os.remove(f"{chrom_fa_path}.fai")
    
    with open(prefix + '.read1.fastq.gz', 'w+') as f1, open(prefix + '.read2.fastq.gz', 'w+') as f2:
        call = subprocess.run(['cat', prefix + '_allele0.read1.fastq.gz', prefix + '_allele1.read1.fastq.gz'], stdout=f1)
        call = subprocess.run(['cat', prefix + '_allele0.read2.fastq.gz', prefix + '_allele1.read2.fastq.gz'], stdout=f2)
    os.remove(prefix + '_allele0.read1.fastq.gz')
    os.remove(prefix + '_allele1.read1.fastq.gz')
    os.remove(prefix + '_allele0.read2.fastq.gz')
    os.remove(prefix + '_allele1.read2.fastq.gz')

def gen_reads(
    tree: Tree,
    params: Optional[Parameters] = None,
    out_dir: Optional[str] = None,
    ref: Optional[str] = None,
    alt_ref: Optional[str] = None,
    coverage: Optional[float] = None,
    bin_size: Optional[int] = None,
    read_len: Optional[int] = None,
    lorenz_y: Optional[float] = None,
    num_processors: Optional[int] = None,
    lorenz_x: float = 0.5
):
    """Given a Tree object returned by GrowthSimulator.simulate_singlecell_lineage, generates synthetic paired-end DNA sequencing reads for each cell. Requires the *dwgsim* and *samtools* binaries be available in the user's :code:`$PATH` variable. Will create two fastq files for each cell located in the 'reads' directory within the provided output directory. See :ref:`Parameters <parameters>` for an explanation of the parameters.

    Args:
        tree (Tree):
        params (Parameters, optional): 
        out_dir (str, optional): 
        ref (str, optional): 
        alt_ref (str, optional): 
        coverage (float, optional): 
        bin_size (int, optional): 
        read_len (int, optional): 
        lorenz_y (float, optional): 
        num_processors (int, optional): 
    """
    params = fill_params(params, out_dir=out_dir, ref=ref, alt_ref=alt_ref, coverage=coverage, bin_size=bin_size, read_len=read_len, lorenz_y=lorenz_y, num_processors=num_processors)

    if params.ref is None:
        raise ValueError(f"Must pass a reference genome.")

    read_dir = os.path.join(params.out_dir, 'reads')
    if not os.path.isdir(read_dir):
        os.makedirs(read_dir)

    check_dependencies(['samtools', 'dwgsim'])

    readcounts = gen_readcounts_singlecell(tree, params=params, lorenz_x=lorenz_x)
    leaves = list(tree.iter_leaves())

    if params.num_processors == 1:
        for cell in leaves:
            gen_reads_cell(cell, readcounts[cell], params.out_dir, params.ref, params.alt_ref, params.read_len, params.seq_error)
    else:
        with multiprocessing.Pool(processes=params.num_processors) as pool:
            args = [(cell, readcounts[cell], params.out_dir, params.ref, params.alt_ref, params.read_len, params.seq_error) for cell in leaves]
            pool.starmap(gen_reads_cell, args)
            
def gen_readcounts_singlecell(
    tree: Tree,
    params: Optional[Parameters] = None,
    out_dir: Optional[str] = None,
    coverage: Optional[float] = None,
    bin_size: Optional[int] = None,
    read_len: Optional[int] = None,
    lorenz_y: Optional[float] = None,
    lorenz_x: float = 0.5
):
    """Given a Tree object returned by GrowthSimulator.simulate_singlecell_lineage, compiles cell-specific region-level copy numbers and draws allele-specific read counts. These will automatically be saved in the output directory. See :ref:`Parameters <parameters>` for an explanation of the parameters.

    Args:
        tree (Tree):
        params (Parameters, optional): 
        out_dir (str, optional): 
        coverage (float, optional): 
        bin_size (int, optional): 
        read_len (int, optional): 
        lorenz_y (float, optional): 

    Returns:
        dict: A dictionary where keys are cells and values are a dictionary containing the region-level allele-specific copy numbers.
    """
    params = fill_params(params, out_dir=out_dir, coverage=coverage, bin_size=bin_size, read_len=read_len, lorenz_y=lorenz_y)
    region_len = tree.root.library.region_len
    regions = tree.root.library.regions
    chrom_lens = tree.root.library.chrom_lens
    mean_rc = (region_len*(params.coverage/2)) / (2*params.read_len)
    regions_per_bin = params.bin_size / region_len
    [Aa, Bb] = get_alpha_beta(lorenz_x, params.lorenz_y)

    last_region_lens = {chrname: int(l % region_len) for chrname,l in chrom_lens.items()}

    leaves = list(tree.iter_leaves())
    gt_readcounts = {}
    region_readcounts = {cell: {} for cell in leaves}

    for cell in leaves:
        full_readcounts = {chrname: [[0 for i in range(regions[chrname])], [0 for i in range(regions[chrname])]] for chrname in chrom_lens}
        for chrname in chrom_lens:
            for chromosome in cell.genome[chrname]:
                if params.lorenz_y == 0.5: # Uniform Coverage
                    cur_readcounts = [round(poisson.rvs(mean_rc)) for i in range(len(chromosome.seq))]

                else: # Non-uniform Coverage
                    cur_readcounts = draw_readcounts_nonuniform(len(len(chromosome.seq)), Aa, Bb, mean_rc, interval=5)
                region_readcounts[cell][chromosome] = cur_readcounts

                for i,r in enumerate(chromosome.seq):
                    q = get_reg_id(r)
                    rc = cur_readcounts[i]
                    if q == regions[chrname] - 1:
                        rc = rc * (last_region_lens[chrname] / region_len)
                    full_readcounts[chrname][chromosome.allele][q] += rc


        binned_readcounts = {chrname: [] for chrname in chrom_lens}
        for chrname in chrom_lens:
            for allele in [0, 1]:
                binned_readcounts[chrname].append(bin_region_readcounts(full_readcounts[chrname][allele], regions_per_bin))
        gt_readcounts[cell] = binned_readcounts
    save_singlecell_readcounts(gt_readcounts, params.out_dir, tree.root.library.chrom_lens, params.bin_size)

    return region_readcounts

def gen_readcounts_bulk(
    tree: Tree,
    params: Optional[Parameters] = None,
    out_dir: Optional[str] = None,
    coverage: Optional[float] = None
):
    """Given a Tree object returned by GrowthSimulator.simulate_clonal_lineage, draws a total and variant read count for each site at every SNV position in the observed cells. These will automatically be saved in the output directory. See :ref:`Parameters <parameters>` for an explanation of the parameters.


    Args:
        tree (Tree):
        params (Parameters, optional): 
        out_dir (str, optional): 
        coverage (float, optional): 
    """
    params = fill_params(params, out_dir=out_dir, coverage=coverage)
    chrom_lens = tree.root.library.chrom_lens
    region_len = tree.root.library.region_len
    chrnames = sort_chrom_names(chrom_lens.keys())

    observed_clones = [node for node in tree.iter_leaves()]
    mutated_basepairs = get_mutated_basepairs(observed_clones, chrnames, region_len)

    nsites = len(set([clone.site for clone in observed_clones]))

    for clone in observed_clones:
        if 'count' not in clone.info:
            raise ValueError("Leaves in Tree have not been initialized with a 'count' key/pair in the info attribute. Ensure you have run GrowthSimulator.simulate_clonal_lineage.")
    
    #Count the number of cells sampled from each clone, separated by site
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
    
    #Compute readcounts
    readcounts = {s: [] for s in range(nsites)}
    for s in range(nsites):
        ncells = sum(list(cell_counts[s].values()))
        cur_clones = list(cell_counts[s].keys())
        if ncells == 0:
            continue
        props = [cell_counts[s][clone] / ncells for clone in cur_clones]
        for chrname,pos in mutated_basepairs:
            r = pos // region_len
            p = pos % region_len

            fraction = 0
            CNs = []
            for clone,prop in zip(cur_clones, props):
                x_ref, x_alt = 0, 0
                for chromosome in clone.genome[chrname]:
                    for q in chromosome.seq:
                        if get_reg_id(q) == r:
                            x_ref += 1
                            if q in chromosome.SNVs:
                                if p in chromosome.SNVs[q]:
                                    x_ref -= 1
                                    x_alt += 1
                fraction += (prop * (x_alt / (x_ref + x_alt + 1e-8)))
                CNs.append(x_ref + x_alt)

            scale = np.mean(CNs) / 2
            exp_rc = params.coverage * scale
            tot_rc = poisson.rvs(exp_rc)
            var_rc = round(tot_rc*fraction)
            readcounts[s].append([tot_rc, var_rc])
    
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

    save_clonal_readcounts(mutated_basepairs, readcounts, site_ids, params.out_dir)


