import numpy as np
import copy
import os
from typing import Optional, Union, Dict

from sistem.selection.base_library import BaseLibrary, Attributes
from sistem.anatomy.utils import assign_n_alt_drivers
from sistem.utilities.utilities import get_reg_id
from sistem.parameters import Parameters, fill_params

class BaseRegionLibrary(BaseLibrary):
    """
    Base region library

    Attributes:
        drivers (dict): keys are chromosome names and values are a list of ids of OGs and TSGs
        essential (dict): keys are chromosome names and values are a list of ids of ESSs
    """
    is_driver_region_model = True
    is_driver_SNV_model = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.drivers = {chrname: [] for chrname in self.regions.keys()}
        self.essential = {chrname: [] for chrname in self.regions.keys()}
    
    # Needed
    def compute_fitness(self, clone):
        """
        Update the fitness coefficient of the genome by traversing over every key region
        """
        fitness = 1
        ploidy = clone.get_ploidy()
        for chrname in self.regions:
            for r in self.drivers[chrname]:
                s = self[chrname][r]
                n = clone.attributes.driver_counts[chrname][r]
                fitness *= (1 + s)**(n/ploidy)
        return fitness
    
    # Needed
    def init_base_fit(self):
        self.base_fit = 1
        for chrname, regions in self.drivers.items():
            for r in regions:
                s = self[chrname][r]
                self.base_fit *= (1 + s)
    
    def init_max_fit(self):
        def compute_fitness_from_counts(driv_counts, ploidy):
            if ploidy > self.max_ploidy or ploidy < self.min_ploidy:
                return -1
            fitness = 1
            for chrname in self.regions:
                for r,n in driv_counts[chrname].items():
                    s = self[chrname][r]
                    fitness *= (1 + s)**(n/ploidy)
            return fitness

        OG, TSG = [], []
        for chrname,deltas in self.delta.items():
            for r,s in enumerate(deltas):
                if s > 0:
                    OG.append((chrname,r,s))
                elif s < 0:
                    TSG.append((chrname,r,s))
        OG.sort(key=lambda x: x[2], reverse=True)
        TSG.sort(key=lambda x: x[2])

        driv_counts = {chrname: {r: 2 for r in self.drivers[chrname]} for chrname in self.regions}
        cur_fit = self.base_fit
        start_nregions = self._ndiploid_regions
        cur_nregions = start_nregions*2
        count = 0
        while count < self.max_distinct_driv and (len(OG) > 0 or len(TSG) > 0):
            fit1, fit2 = -1, -1
            if len(OG) > 0:
                chrname,r,s = OG[0]
                driv_counts[chrname][r] = self.max_region_CN
                nregions = cur_nregions + self.max_region_CN - 2
                fit1 = compute_fitness_from_counts(driv_counts, nregions/start_nregions)
                driv_counts[chrname][r] = 2
            if len(TSG) > 0:
                chrname,r,s = TSG[0]
                driv_counts[chrname][r] = 0
                nregions = cur_nregions - 2
                fit2 = compute_fitness_from_counts(driv_counts, nregions/start_nregions)
                driv_counts[chrname][r] = 2
            if fit1 == -1 and fit2 == -1:
                break
            elif fit1 > fit2:
                chrname,r,s = OG.pop(0)
                driv_counts[chrname][r] = self.max_region_CN
                cur_nregions = cur_nregions + self.max_region_CN - 2
                cur_fit = fit1
            else:
                chrname,r,s = TSG.pop(0)
                driv_counts[chrname][r] = 0
                nregions = cur_nregions - 2
                cur_fit = fit2
            count += 1
        self.max_fit = cur_fit

    def check_viability(self, clone):
        p = 2**(clone.genome.nWGD + 1)
        
        driv_counts = [c for chrname in self.regions for r,c in clone.attributes.driver_counts[chrname].items() if c != p]
        nmutated_drivers = len(driv_counts)
        if len(driv_counts) == 0:
            max_CN = p
        else:
            max_CN = max(driv_counts)
        ploidy = clone.get_ploidy()

        if nmutated_drivers > self.max_distinct_driv:
            return False
        elif ploidy > self.max_ploidy:
            return False
        elif ploidy < self.min_ploidy:
            return False
        elif max_CN > self.max_region_CN:
            return False
        else:
            return True
    
    # Needed
    def update_stats(self, clone, chromosome, start, end, mag=1):
        for r in chromosome.seq[start:end]:
            if r in self.drivers[chromosome.name] + self.essential[chromosome.name]:
                clone.attributes.driver_counts[chromosome.name][r] += mag
    
    # Needed
    def init_attributes(self, clone):
        driver_counts = {chrname: {r: 0 for r in self.drivers[chrname] + self.essential[chrname]} for chrname in self.regions}
        if clone.genome is not None:
            for chrname in self.regions:
                for chromosome in clone.genome[chrname]:
                    for r in chromosome.seq:
                        if r in self.drivers[chrname]:
                            driver_counts[chrname][r] += 1
        clone.attributes = Attributes(driver_counts=driver_counts)
    
    # Needed
    def get_driver_start_regions(self, cell, chromosome, size):
        chrname = chromosome.name
        sequence = chromosome.seq
        possible_start = set()
        for i,r in enumerate(sequence[:-size]):
            # Must be > 0 because it is still in the seq
            if r in self.drivers[chrname] + self.essential[chrname]:
                for j in range(max(0, i - size + 1), i+1):
                    possible_start.add(j)
        return list(possible_start)

    # Needed
    def get_passenger_start_regions(self, cell, chromosome, size):
        """
        Returns a list of all possible start locations such that an event of a given size will not hit any driver regions
        """
        chrname = chromosome.name
        sequence = [get_reg_id(r) for r in chromosome.seq]
        n = len(sequence) - size + 1
        
        starts = []
        i = 0

        while i < n:
            hits = [j for j,r in enumerate(sequence[i:i+size], i) if r in self.drivers[chrname] + self.essential[chrname]]
            if len(hits) > 0:
                i = max(hits) + 1
            else:
                while i < n and sequence[i+size-1] not in self.drivers[chrname] + self.essential[chrname]:
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
        init_drivers = [(chrname, r) for chrname, regions in self.drivers.items() for r in regions]

        ndrivers = len(init_drivers)
        met_libraries = [self]
        for _ in range(1, nsites):
            new_lib = copy.deepcopy(self)
            n_alt_drivers = min(np.random.binomial(ndrivers, alter_prop), ndrivers)
            idxs = np.random.choice(range(ndrivers), n_alt_drivers, replace=False)
            alt_drivers = [init_drivers[i] for i in idxs]
            for chrname, r in alt_drivers:
                new_lib.delta[chrname][r] = np.random.uniform(-CN_coeff, CN_coeff)
            new_lib.init_base_fit()
            new_lib.init_max_fit()
            met_libraries.append(new_lib)

        return met_libraries
    
    def create_metastatic_libraries_random_alt_drivers(self, nsites, alter_prop, CN_coeff, OG_r, TSG_r):
        """
        Similar to create_metastatic_libraries_random, creates metastatic site libraries by randomly altering a N driver regions, where N is drawn from a binomial with a mean of alter_prop. However, instead of resampling the coefficients of the chosen N drivers, it sets the coefficients to 0 and selects N random regions from among all possible.
        """
        OG_p = OG_r / (OG_r + TSG_r)

        prim_drivers = [(chrname, r) for chrname, regions in self.drivers.items() for r in regions]
        prim_neutrals = [(chrname, r) for chrname, num_regions in self.regions.items() for r in range(num_regions) if r not in self.drivers[chrname]]

        num_drivers = len(prim_drivers)
        num_neutrals = len(prim_neutrals)
        met_libraries = [self]
        for _ in range(1, nsites):
            new_lib = copy.deepcopy(self)
            num_former_drivers = np.random.binomial(int(num_drivers*alter_prop), 0.5)
            num_new_drivers = num_drivers - num_former_drivers

            former_driver_idxs = np.random.choice(range(num_drivers), num_former_drivers, replace=False)
            new_driver_idxs = np.random.choice(range(num_neutrals), num_new_drivers, replace=False)

            for i in former_driver_idxs:
                chrname, r = prim_drivers[i]
                new_lib.delta[chrname][r] = 0
                new_lib.drivers[chrname].remove(r)
            
            for i in new_driver_idxs:
                chrname, r = prim_neutrals[i]
                new_lib.drivers[chrname].append(r)
                if np.random.binomial(1, OG_p) == 1:
                    new_lib.delta[chrname][r] = np.random.uniform(0, CN_coeff)
                else:
                    new_lib.delta[chrname][r] = np.random.uniform(-CN_coeff, 0)

            for chrname in new_lib.drivers:
                new_lib.drivers[chrname].sort()

            new_lib.init_base_fit()
            new_lib.init_max_fit()
            met_libraries.append(new_lib)

        return met_libraries
    
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
           
class RandomRegionLibrary(BaseRegionLibrary):
    """The region selection model with randomly generated selection coefficients.
    
    This model is similar to the chromosome-arm model but operates at the region/gene level. In particular, each reference region i on chromosome k is assigned a selection coefficient in range (-CN_coeff, CN_coeff), where positive values correspond to OGs, negative values correspond to tumor suppressor genes TSGs, and zero corresponds to neutral regions (NEU). The fitness :math:`s_a(c)` of cell :math:`c` in site :math:`a` is computed as

    .. math::
    
        s_a(c) = \\prod_{k=1}^K \\prod_{i=1}^{m_k} (1 + \\delta_{k,i})^{x_{k,i} / p_c},

    where :math:`K` is the number of chromosomes, :math:`m_k` is the number of reference regions on chromosome :math:`k`, the :math:`\\delta`'s are selection coefficients, :math:`x_{k,i}` is the copy number of region :math:`i` on chromosome :math:`k` in cell  :math:`c`, and :math:`p_c` is cell ploidy.

    """
    def initialize(
        self, 
        params: Optional[Parameters] = None, 
        CN_coeff: Optional[float] = None,
        OG_r: Optional[float] = None, 
        TSG_r: Optional[float] = None, 
        EG_r: Optional[float] = None
    ):
        """Method used to initialize selection coefficients. See :ref:`Parameters <parameters>` for an explanation of the parameters.

        Args:
            params (Parameters, optional):
            CN_coeff (float, optional):
            OG_r (float, optional):
            TSG_r (float, optional):
            EG_r (float, optional):
        """
        params = fill_params(params, CN_coeff=CN_coeff, OG_r=OG_r, TSG_r=TSG_r, EG_r=EG_r)

        NEU_r = 1 - params.OG_r - params.TSG_r - params.EG_r
        for chrname, num_regions in self.regions.items():
            labels = np.random.choice([0, 1, 2, 3], size=num_regions, p=[params.OG_r, params.TSG_r, params.EG_r, NEU_r])
            for r,l in enumerate(labels):
                if l == 0:
                    self.delta[chrname].append(np.random.uniform(0, params.CN_coeff))
                    self.drivers[chrname].append(r)
                elif l == 1:
                    self.delta[chrname].append(np.random.uniform(-params.CN_coeff, 0))
                    self.drivers[chrname].append(r)
                elif l == 2:
                    self.delta[chrname].append(0)
                    self.essential[chrname].append(r)
                else:
                    self.delta[chrname].append(0)
        self.max_distinct_driv = int(self.max_distinct_driv_ratio * sum([len(self.drivers[chrname]) for chrname in self.drivers]))
        self.init_base_fit()
        self.init_max_fit()
    
class FittedRegionLibrary(BaseRegionLibrary):
    """Uses the same model as :code:`RandomRegionLibrary`, but the selection coefficients are given by the user.

    """
    def initialize(
        self, 
        filepath: Optional[str] = None,
        delta: Optional[Dict] = None,
    ):
        """Method used to initialize selection coefficients.

        Args:
            filepath (str, optional): Path to a tsv file containing the selection coefficients. The first entry is the chromosome name, followed by the selection coefficients in order of reference regions.
            delta (dict, optional): Instead of passing a file, users can pass a dictionary where keys are chromosome names and values are an :math:`m_k`-length array of coefficients.
        
        """
        if delta is not None:
            if len(set(delta.keys()) - set(self.chrom_lens.keys())) > 0:
                raise ValueError("Provided selection coefficient dictionary contains unrecognized chromosomes.")
            
            deltas = {}
            for chrname,v in delta.items():
                if not isinstance(v, (list, tuple, np.ndarray)):
                    raise ValueError("Provided selection coefficient dictionary values must be an array or list containing float values.")
                if len(v) > self.regions[chrname]:
                    raise ValueError(f"Number of deltas for {chrname} exceeds the number of regions specified by chrom_lens and region_len on library initialization.")
                for x in v:
                    if not isinstance(x, (float, np.floating)):
                        raise ValueError(f"Provided deltas are neither floats nor np.floats.")
                v = [float(x) for x in v]
                while len(v) != self.regions[chrname]:
                    v.append(0)
                deltas[chrname] = v
        
        elif filepath is not None:
            assert os.path.isfile(filepath), 'File does not exist at specified path.'
            deltas = {}
            with open(filepath) as f:
                for line in f:
                    info = line.strip().split('\t')
                    chrname = info[0]
                    if chrname not in self.chrom_lens:
                        raise ValueError("Provided arm selection coefficients file contains unrecognized chromosome.")
                    deltas[chrname] = [float(x) for x in info[1:]]

        else:
            raise ValueError("Must pass either filepath or delta.")
        
        for chrname in self.regions:
            if chrname not in deltas:
                deltas[chrname] = [0 for _ in range(self.regions[chrname])]
        self.delta = deltas
        
        for chrname,v in deltas.items():
            for i,x in enumerate(v):
                if x != 0:
                    self.drivers[chrname].append(i)

        self.max_distinct_driv = int(self.max_distinct_driv_ratio * sum([len(self.drivers[chrname]) for chrname in self.drivers]))
        self.init_base_fit()
        self.init_max_fit()