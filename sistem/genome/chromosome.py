import numpy as np
import copy
from collections import defaultdict

from sistem.utilities.utilities import get_reg_id, update

class Chromosome:
    def __init__(self, name='', allele=None, homolog_id=None, seq=None, cent=None):
        self.name = name
        self.allele = allele
        self.homolog_id = homolog_id
        self.seq = seq
        self.cent = cent

    def __len__(self):
        return len(self.seq)
    
    def intact(self):
        return len(self.seq)

    def copy(self, homolog_id=None):
        if homolog_id == None:
            homolog_id = self.homolog_id
        new_chrom = Chromosome(name=self.name, allele=self.allele, homolog_id=homolog_id, seq=copy.copy(self.seq), cent=self.cent)
        return new_chrom
    
    # place m additional copies of the segment from start to end in tandem after m
    def amplify(self, start, end, m):
        new_seg = m * self.seq[start:end]
        self.seq = self.seq[:end] + new_seg + self.seq[end:]
        if self.cent is not None and self.cent >= end:
            self.cent += m*(end - start)
    
    # Delete the segment from start to end
    def delete(self, start, end):
        self.seq = self.seq[:start] + self.seq[end:]
        if self.cent is not None:
            if self.cent > start and self.cent < end:
                self.cent = start
            if self.cent >= end:
                self.cent -= (end - start)
            if self.cent == 0 or self.cent >= len(self.seq):
                self.cent = None



def default_SNV_obj():
    return defaultdict(int)

class SNVChromosome(Chromosome):
    """
    Attributes:
        unique_regions (defaultdict(int)): counts the number of times some region r has appeared to create unique region id
        SNVs (defaultdict(defaultdict(int))): dictionary of dictionary of ints. First key is region id; Second key is the position in that region (between 0 and region_len - 1), third is either 0,1,2, referring to the non-reference base pair.

    """
    def __init__(self, unique_regions=None, SNVs=None, driverSNV_counts=None, ndriverSNVs=None, **kwargs):
        self.unique_regions = unique_regions if unique_regions is not None else defaultdict(int)
        self.SNVs = SNVs if SNVs is not None else defaultdict(default_SNV_obj)
        #self.driverSNV_counts = driverSNV_counts if driverSNV_counts is not None else defaultdict(int)
        super().__init__(**kwargs)

    def copy(self, homolog_id=None):
        if homolog_id == None:
            homolog_id = self.homolog_id
        new_chrom = SNVChromosome(name=self.name, allele=self.allele, homolog_id=homolog_id, seq=copy.copy(self.seq), cent=self.cent, unique_regions=copy.copy(self.unique_regions), SNVs=copy.copy(self.SNVs))#, driverSNV_counts = copy.copy(self.driverSNV_counts))
        return new_chrom
    
    def up(self, r):
        v = get_reg_id(r)
        self.unique_regions[v] += 1
        return update(v, self.unique_regions[v])

    def amplify(self, start, end, m):
        gain_seg = self.seq[start:end]
        #driverSNV_regions = [i for i,r in enumerate(gain_seg) if r in self.driverSNV_counts]
        SNV_regions = [i for i,r in enumerate(gain_seg) if r in self.SNVs]
        new_seg = []
        for _ in range(m):
            temp_new_seg = [self.up(r) for r in gain_seg]
            new_seg.extend(temp_new_seg)
            #for i in driverSNV_regions:
            #    self.driverSNV_counts[temp_new_seg[i]] = np.copy(self.driverSNV_counts[gain_seg[i]])
            for i in SNV_regions:
                self.SNVs[temp_new_seg[i]] = copy.copy(self.SNVs[gain_seg[i]])
        self.seq = self.seq[:end] + new_seg + self.seq[end:]
        if self.cent is not None and self.cent >= end:
            self.cent += m*(end - start)
    
    def delete(self, start, end):
        for r in self.seq[start:end]:
            #if r in self.driverSNV_counts:
            #    del self.driverSNV_counts[r]
            if r in self.SNVs:
                del self.SNVs[r]
        super().delete(start, end)
    
    def add_SNV(self, start, position, bp):
        self.SNVs[int(self.seq[start])][position] = bp