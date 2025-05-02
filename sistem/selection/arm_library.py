import numpy as np
import copy
import os
from collections import Counter
from typing import Optional, Union, Dict

from sistem.selection.base_library import BaseLibrary, Attributes
from sistem.anatomy.utils import assign_n_alt_drivers
from sistem.utilities.utilities import get_reg_id
from sistem.parameters  import Parameters, fill_params

class BaseArmLibrary(BaseLibrary):
    is_driver_region_model = False
    is_driver_SNV_model = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._arm_sizes = {}
        for chrname in self.regions.keys():
            if isinstance(self.arm_ratios, dict):
                arm_ratio = self.arm_ratios[chrname]
            else:
                arm_ratio = self.arm_ratios
            cent_idx = round(self.regions[chrname]*arm_ratio)
            # 'p' entry in _arm_sizes is also the cent index
            self._arm_sizes[chrname] = [cent_idx, self.regions[chrname] - cent_idx]
    
    def compute_fitness(self, clone):
        fitness = 1
        ploidy = clone.get_ploidy()
        for chrname in self.regions:
            for a in [0, 1]:
                s = self[chrname][a]
                n = clone.attributes.arm_counts[chrname][a] / self._arm_sizes[chrname][a]
                fitness *= (1 + s)**(n/ploidy)
        return fitness
    
    def init_base_fit(self):
        self.base_fit = 1
        for chrname in self.regions:
            for a in [0, 1]:
                s = self[chrname][a]
                self.base_fit *= (1 + s)

    def init_max_fit(self):
        def compute_fitness_from_counts(arm_counts, ploidy):
            if ploidy > self.max_ploidy or ploidy < self.min_ploidy:
                return -1
            fitness = 1
            for chrname in self.regions:
                for a in [0, 1]:
                    s = self[chrname][a]
                    n = arm_counts[chrname][a] / self._arm_sizes[chrname][a]
                    fitness *= (1 + s)**(n/ploidy)
            return fitness

        pos, neg = [], []
        for chrname,deltas in self.delta.items():
            for a,s in enumerate(deltas):
                if s > 0:
                    pos.append((chrname,a,s))
                elif s < 0:
                    neg.append((chrname,a,s))
        pos.sort(key=lambda x: x[2], reverse=True)
        neg.sort(key=lambda x: x[2])

        arm_counts = {chrname: [2*s for s in sizes] for chrname, sizes in self._arm_sizes.items()}
        cur_fit = self.base_fit
        start_nregions = self._ndiploid_regions
        cur_nregions = start_nregions*2
        count = 0
        while count < self.max_distinct_driv and (len(pos) > 0 or len(neg) > 0):
            fit1, fit2 = -1, -1
            if len(pos) > 0:
                chrname,a,s = pos[0]
                arm_counts[chrname][a] = self.max_region_CN*self._arm_sizes[chrname][a]
                nregions = cur_nregions + arm_counts[chrname][a] - (2*self._arm_sizes[chrname][a])
                fit1 = compute_fitness_from_counts(arm_counts, nregions/start_nregions)
                arm_counts[chrname][a] = 2*self._arm_sizes[chrname][a]
            if len(neg) > 0:
                chrname,a,s = neg[0]
                arm_counts[chrname][a] = 0
                nregions = cur_nregions - (2*self._arm_sizes[chrname][a])
                fit2 = compute_fitness_from_counts(arm_counts, nregions/start_nregions)
                arm_counts[chrname][a] = 2*self._arm_sizes[chrname][a]
            if fit1 == -1 and fit2 == -1:
                break
            elif fit1 > fit2:
                chrname,a,s = pos.pop(0)
                arm_counts[chrname][a] = self.max_region_CN*self._arm_sizes[chrname][a]
                cur_nregions = cur_nregions + arm_counts[chrname][a] - (2*self._arm_sizes[chrname][a])
                cur_fit = fit1
            else:
                chrname,a,s = neg.pop(0)
                arm_counts[chrname][a] = 0
                cur_nregions = cur_nregions - (2*self._arm_sizes[chrname][a])
                cur_fit = fit2
            count += 1
        self.max_fit = cur_fit
    
    def check_viability(self, clone):
        """
        Checks if clone passes viability checkpoints based on mutated driver stats. 

        Args:
            clone (Clone): The clone object in question.

        Returns:
            (bool): True if passes, False if not.
        
        """
        max_chrom_CNs = [0]
        for chrname in self.regions:
            c = Counter()
            for chromosome in clone.genome[chrname]:
                c.update(Counter(chromosome.seq))
            if len(c):
                max_chrom_CNs.append(c.most_common()[0][1])
        max_CN = max(max_chrom_CNs)
        
        p = 2**(clone.genome.nWGD + 1)
        nmutated_arms = 0
        for chrname in self.regions:
            for a in [0, 1]:
                r = clone.attributes.arm_counts[chrname][a]/(self._arm_sizes[chrname][a]*p)
                if r > 1.5 or r < 0.5:
                    nmutated_arms += 1

        ploidy = clone.get_ploidy()

        if nmutated_arms > self.max_distinct_driv:
            return False
        elif ploidy > self.max_ploidy:
            return False
        elif ploidy < self.min_ploidy:
            return False
        elif max_CN > self.max_region_CN:
            return False
        else:
            return True

    def update_stats(self, clone, chromosome, start, end, mag=1):
        for r in chromosome.seq[start:end]:
            if r < self._arm_sizes[chromosome.name][0]:
                clone.attributes.arm_counts[chromosome.name][0] += mag
            else:
                clone.attributes.arm_counts[chromosome.name][1] += mag
    
    def init_attributes(self, clone):
        arm_counts = {chrname: [0, 0] for chrname in self.regions}
        if clone.genome is not None:
            for chromosome in clone.genome.get_chromosomes():
                chrname = chromosome.name
                for r in chromosome.seq:
                    if r < self._arm_sizes[chrname][0]:
                        arm_counts[chrname][0] += 1
                    else:
                        arm_counts[chrname][1] += 1
        clone.attributes = Attributes(arm_counts=arm_counts)

    def get_driver_start_regions(self, cell, chromosome, size):
        n = len(chromosome) - size
        return list(range(n))

    def get_passenger_start_regions(self, cell, chromosome, size):
        sequence = [get_reg_id(r) for r in chromosome.seq]
        n = len(sequence) - size + 1
        blregions = None
        if chromosome.name in cell.attributes.blacklisted_regions:
            if chromosome.homolog_id in cell.attributes.blacklisted_regions[chromosome.name]:
                blregions = set(cell.attributes.blacklisted_regions[chromosome.name][chromosome.homolog_id])
        if blregions is None:
            return list(range(n))

        starts = []
        i = 0

        while i < n:
            hits = [j for j,r in enumerate(sequence[i:i+size], i) if r in blregions]
            if len(hits) > 0:
                i = max(hits) + 1
            else:
                while i < n and sequence[i+size-1] not in blregions:
                    starts.append(i)
                    i += 1
                    
        return starts
    
    def create_metastatic_libraries(self, nsites, alter_prop, CN_coeff, method='random', dists=None):
        if method == 'random':
            return self.create_metastatic_libraries_random(nsites, alter_prop, CN_coeff)
        elif method == 'distance':
            if dists is None:
                raise ValueError("Cannot create metastatic libraries with the 'distance' method when dists has not been initialized.")
            return self.create_metastatic_libraries_distance(nsites, alter_prop, CN_coeff, dists)
        else:
            raise ValueError("Method not recognized.")
    
    def create_metastatic_libraries_random(self, nsites, alter_prop, CN_coeff):
        init_drivers = [(chrname, a) for chrname in self.regions for a in [0, 1]]

        ndrivers = len(init_drivers)
        libraries = [self]
        for _ in range(1, nsites):
            new_lib = copy.deepcopy(self)
            n_alt_drivers = min(np.random.binomial(ndrivers, alter_prop), ndrivers)
            idxs = np.random.choice(range(ndrivers), n_alt_drivers, replace=False)
            alt_drivers = [init_drivers[i] for i in idxs]
            for chrname, a in alt_drivers:
                new_lib.delta[chrname][a] = np.random.uniform(-CN_coeff, CN_coeff)
            new_lib.init_base_fit()
            new_lib.init_max_fit()
            libraries.append(new_lib)

        return libraries

    def create_metastatic_libraries_distance(self, nsites, alter_prop, CN_coeff, dists):
        init_drivers = [(chrname, a) for chrname in self.regions for a in [0, 1]]

        ndrivers = len(init_drivers)
        n_alt_drivers = np.random.binomial(ndrivers, alter_prop)
        dists = np.rint(dists*n_alt_drivers).astype(int)
        idxs = np.triu_indices(nsites)
        dmatrix = np.zeros((nsites, nsites), dtype=int)
        dmatrix[idxs] = dists[:, 0]
        dmatrix = dmatrix + dmatrix.T
        A, score = assign_n_alt_drivers(dmatrix, ndrivers)

        libraries = [self]
        unique = list(set.union(*[x for x in A.values()]))
        idxs = np.random.choice(range(ndrivers), len(unique), replace=False)
        alt_drivers = [init_drivers[i] for i in idxs]
        new_drivers = {}
        for i,(chrname, a) in zip(unique, alt_drivers):
            new_delta = np.random.uniform(-CN_coeff, CN_coeff)
            new_drivers[i] = (chrname, a, new_delta)

        for s in range(1, nsites):
            new_lib = copy.deepcopy(self)
            for i in A[s]:
                chrname, a, delta = new_drivers[i]
                new_lib.delta[chrname][a] = delta
            new_lib.init_base_fit()
            new_lib.init_max_fit()
            libraries.append(new_lib)

        return libraries

class RandomArmLibrary(BaseArmLibrary):
    """The chromosome-arm selection model with randomly generated selection coefficients.
    
    Under this model, each chromosome-arm is assigned a random selection coefficient in the range (-CN_coeff, CN_coeff). Arms with a positive coefficient correspond to those with a greater balance of OGs to TSGs, while arms with a negative coefficient will have a greater balance of TSGs to OGs. The fitness :math:`s_a(c)` of cell :math:`c` in site :math:`a` is computed as

    .. math::
    
        s_a(c) = \\prod_{k=1}^K \\prod_{arm \\in \\{p,q\\}} (1 + \\delta_{k,arm})^{x_{k,arm} / p_c},

    where :math:`K` is the number of chromosomes, :math:`\\delta` are selection coefficients, :math:`x_{k,arm}` is the average copy number of chromosome arm :math:`k,arm` in cell :math:`c`, and :math:`p_c` is cell ploidy.

    """
    def initialize(
        self, 
        params: Optional[Parameters] = None, 
        CN_coeff: Optional[float] = None
    ):
        """Method used to initialize selection coefficients. See :ref:`Parameters <parameters>` for an explanation of the parameters.

        Args:
            params (Parameters, optional):
            CN_coeff (float, optional):
        
        """
        params = fill_params(params, CN_coeff=CN_coeff)

        arm_deltas = {}
        for chrname in self.regions:
            arm_deltas[chrname] = list(np.random.uniform(-params.CN_coeff, params.CN_coeff, 2))
        self.delta = arm_deltas
        self.max_distinct_driv = int(self.max_distinct_driv_ratio * (2*len(self.regions)))
        self.init_base_fit()
        self.init_max_fit()

class FittedArmLibrary(BaseArmLibrary):
    """Uses the same model as :code:`RandomArmLibrary`, but the selection coefficients are given by the user.
    
    """
    def initialize(
        self, 
        filepath: Optional[str] = None,
        delta: Optional[Dict] = None,
    ):
        """Method used to initialize selection coefficients.

        Args:
            filepath (str, optional): Path to a tsv file containing the selection coefficients. The first entry is the chromosome name, followed by the selection coefficient of the short and long arm, respectively.
            delta (dict, optional): Instead of passing a file, users can pass a dictionary where keys are chromosome names and values are a 2-length array of coefficients.
        
        """
        if delta is not None:
            if len(set(delta.keys()) - set(self.chrom_lens.keys())) > 0:
                raise ValueError("Provided selection coefficient dictionary contains unrecognized chromosomes.")
            for i,v in delta.items():
                if not isinstance(v, (list, tuple, np.ndarray)):
                    raise ValueError("Provided selection coefficient dictionary values must be an array, list, or tuple containing 2 float values.")
                v = list(v)
                if len(v) != 2 or not (isinstance(v[0], (float, np.floating)) or isinstance(v[1], (float, np.floating))):
                    raise ValueError("Provided selection coefficient dictionary values must be an array, list, or tuple containing 2 float values.")
            arm_deltas = {i: list(v) for i,v in delta.items()}
        elif filepath is not None:
            assert os.path.isfile(filepath), 'File does not exist at specified path'
            arm_deltas = {}
            with open(filepath) as f:
                for line in f:
                    info = line.strip().split('\t')
                    chrname = info[0]
                    if chrname not in self.chrom_lens:
                        raise ValueError("Provided arm selection coefficients file contains unrecognized chromosome.")
                    arm_deltas[chrname] = [float(info[1]), float(info[2])]
        else:
            raise ValueError("Must pass either filepath or delta.")
        
        for chrname in self.chrom_lens:
            if chrname not in arm_deltas:
                arm_deltas[chrname] = [0, 0]
        self.delta = arm_deltas
        self.max_distinct_driv = int(self.max_distinct_driv_ratio * (2*len(self.regions)))
        self.init_base_fit()
        self.init_max_fit()
