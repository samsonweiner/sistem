import numpy as np
from collections.abc import Mapping

from sistem.genome.chromosome import Chromosome, SNVChromosome

def init_diploid_genome(library):
    genome = {chrname: [] for chrname in library.regions.keys()}
    for chrname,num_regions in library.regions.items():
        if isinstance(library.arm_ratios, dict):
            arm_ratio = library.arm_ratios[chrname]
        else:
            arm_ratio = library.arm_ratios

        cent_idx = round(num_regions*arm_ratio)
        hap_profile = [i for i in range(num_regions)]
        # If SNVs are involved
        if library.is_driver_SNV_model:
            genome[chrname].append(SNVChromosome(name=chrname, allele=0, homolog_id=0, seq=[i for i in hap_profile], cent = cent_idx))
            genome[chrname].append(SNVChromosome(name=chrname, allele=1, homolog_id=1, seq=[i for i in hap_profile], cent = cent_idx))
        # If SNVs are not involved
        else:
            genome[chrname].append(Chromosome(name=chrname, allele=0, homolog_id=0, seq=[i for i in hap_profile], cent = cent_idx))
            genome[chrname].append(Chromosome(name=chrname, allele=1, homolog_id=1, seq=[i for i in hap_profile], cent = cent_idx))
    g = Genome(genome)
    return g

class Genome(Mapping):
    def __init__(self, genome, nregions=None):
        self.genome = genome
        self.nWGD = 0

        self.nregions = nregions
        if self.nregions == None:
            self.nregions = len(self)

    def __getitem__(self, chrname):
        '''
        Returns the two allele lists of the chromosome.
        '''
        return self.genome[chrname]
    
    def __iter__(self):
        '''
        Returns an iterable over chromosome names.
        '''
        return iter(self.genome)
    
    def __len__(self):
        l = 0
        for chrname in self:
            for chromosome in self[chrname]:
                l += len(chromosome)
        return l
    
    @property
    def chrom_names(self):
        return [l for l in self]
    
    def find(self, chrname, homolog_id):
        for chromosome in self[chrname]:
            if chromosome.homolog_id == homolog_id:
                return chromosome
    
    def get_chromosomes(self):
        flat_chroms = [chromosome for chrname in self for chromosome in self[chrname]]
        return flat_chroms

    def get_random_chrom(self, exclude=[]):
        '''
        Returns a random chromosome with probability proportional to its length, excluding those passed by argument.

        Args:
            exclude (list): List of Chromosome instances to exclude from selection
        
        Returns:
            chromosome (Optional[Chromosome, None]): A randomly selected chromosome if one is available, otherwise None.
        '''
        flat_chroms = [chromosome for chrname in self for chromosome in self[chrname] if chromosome not in exclude]
        chrom_sizes = [len(chromosome) for chromosome in flat_chroms]
        l = sum(chrom_sizes)
        if len(flat_chroms) == 0 or l == 0:
            return None
        p = [len(chromosome) / l for chromosome in flat_chroms]
        chromosome = np.random.choice(flat_chroms, p=p)
        return chromosome

    def get_random_intact_chrom(self):
        '''
        Returns a random chromosome which still has both arms intact
        '''
        intact_chroms = []
        for chrname in self:
            for chromosome in self[chrname]:
                if chromosome.intact() and chromosome.cent:
                    intact_chroms.append(chromosome)
        chromosome = np.random.choice(intact_chroms)
        return chromosome

    def get_random_intact_arm(self):
        '''
        Returns a chromosome arm from among all fully intact or partially intact chromosomes
        '''
        intact_arms = []
        for chrname in self:
            for chromosome in self[chrname]:
                if chromosome.intact():
                    if chromosome.cent:
                        intact_arms.append((chromosome, 'p'))
                        intact_arms.append((chromosome, 'q'))
                    else:
                        intact_arms.append((chromosome, 'w'))
        if len(intact_arms) > 0:
            idx = np.random.choice(list(range(len(intact_arms))))
            chromosome, arm = intact_arms[idx]
            return (chromosome, arm)
        else:
            return None