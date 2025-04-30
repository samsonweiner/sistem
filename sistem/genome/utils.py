from pyfaidx import Fasta
import copy

from sistem.genome.chromosome import SNVChromosome
from sistem.genome.genome import Genome

def ConvertToSNVGenome(igenome):
    '''
    Creates an identical genome made of SNVChromosome objects from an initial genome made of standard Chromosome objects

    Args:
        igenome (Genome): Starting genome made of Chromosome objects.

    Returns:
        altgenome (Genome): New genome made of SNVChromosome objects.
    '''
    altgenome = {chrname: [] for chrname in igenome}
    for chrname in altgenome:
        for chromosome in igenome[chrname]:
            altgenome[chrname].append(SNVChromosome(name=chrname, allele=chromosome.allele, homolog_id=chromosome.homolog_id, cent=chromosome.cent, seq=copy.copy(chromosome.seq)))
    g = Genome(altgenome)
    return g

# Gets the number of regions per chrom
def get_num_regions(chrom_lens, region_len):
    regions = {}
    for chrname in chrom_lens:
        nregions = int(chrom_lens[chrname] // region_len)
        rem = nregions % region_len
        if rem > 0:
            nregions += 1
        regions[chrname] = nregions
    return regions

def fixed_chrom_lens(nchroms, chrom_len):
    return {f'chr{i}': chrom_len for i in range(1, nchroms + 1)}

def hg38_chrom_lengths_from_cytoband(include_allosomes=False):
    chrom_lens = {'chr1': 249250621, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr2': 243199373, 'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chrX': 155270560, 'chrY': 59373566}
    arm_ratios = {'chr1': 0.5015, 'chr10': 0.2966, 'chr11': 0.39776, 'chr12': 0.26746, 'chr13': 0.15542, 'chr14': 0.16395, 'chr15': 0.18531, 'chr16': 0.40507, 'chr17': 0.29558, 'chr18': 0.22029, 'chr19': 0.44817, 'chr2': 0.38364, 'chr20': 0.43633, 'chr21': 0.27426, 'chr22': 0.28652, 'chr3': 0.45954, 'chr4': 0.26366, 'chr5': 0.26753, 'chr6': 0.35649, 'chr7': 0.3764, 'chr8': 0.31155, 'chr9': 0.34699, 'chrX': 0.39029, 'chrY': 0.21053}
    if not include_allosomes:
        del chrom_lens['chrX']
        del chrom_lens['chrY']
        del arm_ratios['chrX']
        del arm_ratios['chrY']
    return chrom_lens, arm_ratios

def get_chrom_lens_from_reference(input_fasta, chrom_names=None):
    ref = Fasta(input_fasta, one_based_attributes=False)
    chrom_lens = {}
    for chrom in ref.keys():
        if chrom_names is None:
            chrom_lens[chrom] = len(ref[chrom])
        else:
            if chrom in chrom_names:
                chrom_lens[chrom] = len(ref[chrom])
    return chrom_lens
    