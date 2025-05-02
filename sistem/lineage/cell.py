import copy
import numpy as np
from typing import Optional

from sistem.genome.genome import Genome
from sistem.selection import BaseLibrary, Attributes
from sistem.lineage.mutation import CNA, SNV
from sistem.lineage.tree import Node

class Cell(Node):
    gen = None

    def __init__(
        self, 
        birth_gen: Optional[int] = None, 
        genome: Optional[Genome] = None, 
        library: Optional[BaseLibrary] = None, 
        attributes: Optional[Attributes] = None, 
        site: int = 0, 
        **kwargs
    ):
        self.birth_gen = birth_gen
        self.genome = genome
        self.library = library
        self.attributes = attributes
        self.site = site
        self.events = []
        self.birth_counts = {}
        super().__init__(**kwargs)
    
    def inherit(self):
        assert self.parent is not None, 'Cannot inherit if parent is None'
        self.genome = copy.deepcopy(self.parent.genome)

    def add_SNV(self, chromosome, start, position, bp, record_event=True, driver=False):
        chromosome.add_SNV(start, position, bp)
        if record_event:
            self.events.append(SNV(chrom=chromosome.name, allele=chromosome.allele, homolog_id=chromosome.homolog_id, index=start, region=chromosome.seq[start], position=position, bp=bp, driver=driver, gen=self.gen))

    ## For passengers
    def focal_amplification(self, chromosome, start, end, m, record_event=True, driver=False):
        chromosome.amplify(start, end, m)
        if record_event:
            self.events.append(CNA(category='focal_amplification', chrom=chromosome.name, allele=chromosome.allele, homolog_id=chromosome.homolog_id, start=start, end=end, copies=m, gene_ids=chromosome.seq[start:end], driver=driver, gen=self.gen))
    
    ## For passengers
    def focal_deletion(self, chromosome, start, end, record_event=True, driver=False):
        chromosome.delete(start, end)
        if record_event:
            self.events.append(CNA(category='focal_deletion', chrom=chromosome.name, allele=chromosome.allele, homolog_id=chromosome.homolog_id, start=start, end=end, gene_ids=chromosome.seq[start:end], driver=driver, gen=self.gen))
    
    def arm_duplication(self, chromosome, arm, record_event=True):
        chrname = chromosome.name
        new_homolog_id = len(self.genome[chrname])
        nchromosome = chromosome.copy(homolog_id=new_homolog_id)
        self.genome[chrname].append(nchromosome)
        if nchromosome.cent:
            if arm == 'p':
                nchromosome.delete(nchromosome.cent, len(nchromosome))
            elif arm == 'q':
                nchromosome.delete(0, nchromosome.cent)
        if record_event:
            self.events.append(CNA(category='arm_duplication', chrom=chromosome.name, allele=chromosome.allele, homolog_id=chromosome.homolog_id, spawn_id=new_homolog_id, arm=arm, gen=self.gen))
        
    def arm_deletion(self, chromosome, arm, record_event=True):
        if chromosome.cent:
            if arm == 'p':
                self.focal_deletion(chromosome, 0, chromosome.cent, record_event=False)
            elif arm == 'q':
                self.focal_deletion(chromosome, chromosome.cent, len(chromosome), record_event=False)
        else:
            self.focal_deletion(chromosome, 0, len(chromosome), record_event=False)
        if record_event:
            self.events.append(CNA(category='arm_deletion', chrom=chromosome.name, allele=chromosome.allele, homolog_id=chromosome.homolog_id, arm=arm, gen=self.gen))

    def chromosomal_duplication(self, chromosome, record_event=True):
        chrname = chromosome.name
        new_homolog_id = len(self.genome[chrname])
        self.genome[chrname].append(chromosome.copy(homolog_id=new_homolog_id))
        if record_event:
            self.events.append(CNA(category='chromosomal_duplication', chrom=chromosome.name, allele=chromosome.allele, homolog_id=chromosome.homolog_id, spawn_id=new_homolog_id, gen=self.gen))

    def chromosomal_deletion(self, chromosome, record_event=True):
        self.focal_deletion(chromosome, 0, len(chromosome), record_event=False)
        if record_event:
            self.events.append(CNA(category='chromosomal_deletion', chrom=chromosome.name, allele=chromosome.allele, homolog_id=chromosome.homolog_id, gen=self.gen))

    def WGD(self, record_event=True):
        for chrname in self.genome.chrom_names:
            n = len(self.genome[chrname])
            new_chroms = []
            for chromosome in self.genome[chrname]:
                if len(chromosome) > 0:
                    new_chroms.append(chromosome.copy(homolog_id=n))
                    n += 1
            for nchrom in new_chroms:
                self.genome[chrname].append(nchrom)
        if record_event:
            self.events.append(CNA(category='WGD', gen=self.gen))

    def get_driver_start_regions(self, chromosome, size):
        '''
        Returns a list of all possible start locations such that an event of a given size will hit atleast one driver regions.
        '''
        return self.library.get_driver_start_regions(self, chromosome, size)
    
    def get_passenger_start_regions(self, chromosome, size):
        '''
        Returns a list of all possible start locations such that an event of a given size will not hit any driver regions.
        '''
        return self.library.get_passenger_start_regions(self, chromosome, size)

class Clone(Cell):
    def __init__(
        self, 
        popsize: int = 1, 
        fitness: Optional[float] = None, 
        **kwargs
    ):
        super().__init__(**kwargs)
        self.popsize = popsize
        self.death_gen = None
        self.fitness = fitness
        if self.library and self.attributes is None:
            self.library.init_attributes(self)

    def __len__(self):
        return self.popsize

    def set_library(self, library):
        self.library = library
        self.library.init_attributes(self)
    
    def replicate(self, site=None, popsize=1, library=None, name=None):
        if site == None:
            site = self.site
        if library == None:
            library = self.library

        if library != self.library:
            attributes = self.library.init_attributes(self)
        else:
            attributes = copy.deepcopy(self.attributes)

        new_clone = Clone(name=name, parent=self, genome=copy.deepcopy(self.genome), library=library, site=site, popsize=popsize, attributes=attributes, birth_gen=self.gen, fitness=self.fitness)

        return new_clone
    
    def die(self):
        self.popsize = 0
        self.death_gen = Clone.gen
        # Do not remove genome if object is the starting clone
        if self.birth_gen == 0:
            self.genome = None
            self.attributes = None

    def get_ploidy(self):
        return self.genome.nregions / self.library._ndiploid_regions

    def is_viable(self):
        return self.library.check_viability(self)
    
    def update_fitness(self):
        '''
        Call the library-specific fitness function model
        '''
        self.fitness = self.library.compute_fitness(self)

    def get_birth_prob(self, mean_fit, Ntot, Etot):
        return (self.fitness / mean_fit) * (Etot / (Etot + Ntot))
    
    def get_num_births(self, mean_fit, Ntot, Etot):
        birth_prob = self.get_birth_prob(mean_fit, Ntot, Etot)
        if birth_prob >= 1:
            nbirths = self.popsize
        else:
            nbirths = np.random.binomial(self.popsize, birth_prob)
        return nbirths
    
    def add_SNV(self, chromosome, start, position, bp, record_event=True):
        self.library.update_stats(self, chromosome, start, start)
        super().add_SNV(chromosome, start, position, bp, record_event=record_event, driver=True)

    def focal_amplification(self, chromosome, start, end, m, record_event=True):
        '''
        Given a chromosome, start and end indicies, and a magnitude, amplify the region and update the cell fitness to reflect changes in driver regions
        '''
        self.genome.nregions += m * (end - start)
        self.library.update_stats(self, chromosome, start, end, mag=m)
        super().focal_amplification(chromosome, start, end, m, record_event=record_event, driver=True)

    def focal_deletion(self, chromosome, start, end, record_event=True):
        self.genome.nregions -= (end - start)
        self.library.update_stats(self, chromosome, start, end, mag=-1)
        super().focal_deletion(chromosome, start, end, record_event=record_event, driver=True)

    def arm_duplication(self, chromosome, arm, record_event=True):
        if chromosome.cent:
            if arm == 'p':
                self.genome.nregions += chromosome.cent
                self.library.update_stats(self, chromosome, 0, chromosome.cent)
            elif arm == 'q':
                self.genome.nregions += (len(chromosome) - chromosome.cent)
                self.library.update_stats(self, chromosome, chromosome.cent, len(chromosome))
        else:
            self.genome.nregions += len(chromosome)
            self.library.update_stats(self, chromosome, 0, len(chromosome))
        super().arm_duplication(chromosome, arm, record_event=record_event)
        
    def chromosomal_duplication(self, chromosome, record_event=True):
        self.genome.nregions += len(chromosome)
        self.library.update_stats(self, chromosome, 0, len(chromosome))
        super().chromosomal_duplication(chromosome, record_event=record_event)

    def WGD(self, record_event=True):
        self.genome.nWGD += 1
        self.genome.nregions *= 2
        for chrname in self.genome.chrom_names:
            for chromosome in self.genome[chrname]:
                self.library.update_stats(self, chromosome, 0, len(chromosome))
        super().WGD(record_event=record_event)