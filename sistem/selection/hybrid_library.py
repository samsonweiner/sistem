import numpy as np
import os
from typing import Optional, Dict

from sistem.selection.base_library import Attributes
from sistem.selection.region_library import BaseRegionLibrary
from sistem.utilities.utilities import get_reg_id
from sistem.parameters import Parameters, fill_params

class BaseHybridLibrary(BaseRegionLibrary):
    is_driver_region_model = True
    is_driver_SNV_model = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.lam = {chrname: [] for chrname in self.regions.keys()}
        
    # Needed
    def compute_fitness(self, clone):
        fitness = 1
        ploidy = clone.get_ploidy()
        for chrname in self.regions:
            for r,n in clone.attributes.driver_counts[chrname].items():
                if n == 0:
                    y = 0
                else:
                    y = clone.attributes.SNV_counts[chrname][r]/n
                k = np.exp(-y * self.lam[chrname][r])
                if r in self.drivers[chrname]:
                    s = self[chrname][r]
                    fitness *= (1 + (s*k))**(n/ploidy)
                elif r in self.essential[chrname]:
                    fitness *= k
        return fitness
    
    def check_viability(self, clone):
        p = 2**(clone.genome.nWGD + 1)
        driv_counts, SNV_counts = [], []
        for chrname in self.regions:
            for r,c in clone.attributes.driver_counts[chrname].items():
                if c != p:
                    driv_counts.append(c)
                if c != 0:
                    SNV_counts.append(clone.attributes.SNV_counts[chrname][r]/c)

        nmutated_drivers = len(driv_counts)
        if len(driv_counts) == 0:
            max_CN = p
        else:
            max_CN = max(driv_counts)
        if len(SNV_counts) == 0:
            max_SNV = 0
        else:
            max_SNV = max(SNV_counts)

        ploidy = clone.get_ploidy()

        if nmutated_drivers > self.max_distinct_driv:
            return False
        elif ploidy > self.max_ploidy:
            return False
        elif ploidy < self.min_ploidy:
            return False
        elif max_CN > self.max_region_CN:
            return False
        elif max_SNV >= self.max_region_SNV:
            return False
        else:
            return True

    def update_stats(self, clone, chromosome, start, end, mag=1):
        if start == end:
            v = get_reg_id(start)
            if v in self.drivers[chromosome.name] or v in self.essential[chromosome.name]:
                clone.attributes.SNV_counts[chromosome.name][v] += mag
        for r in chromosome.seq:
            v = get_reg_id(r)
            if v in self.drivers[chromosome.name] or v in self.essential[chromosome.name]:
                clone.attributes.driver_counts[chromosome.name][v] += mag
                clone.attributes.SNV_counts[chromosome.name][v] += mag*len(chromosome.SNVs[r])

    # Needed
    def init_attributes(self, clone):
        driver_counts = {chrname: {r: 0 for r in self.drivers[chrname] + self.essential[chrname]} for chrname in self.regions}
        SNV_counts = {chrname: {r: 0 for r in self.drivers[chrname] + self.essential[chrname]} for chrname in self.regions}
        if clone.genome is not None:
            for chrname in self.regions:
                for chromosome in clone.genome[chrname]:
                    for r in chromosome.seq:
                        v = get_reg_id(r)
                        if v in self.drivers[chrname] + self.essential[chrname]:
                            driver_counts[chrname][v] += 1
                            if r in chromosome.SNVs:
                                SNV_counts[chrname][v] += len(chromosome.SNVs[r])
        clone.attributes = Attributes(driver_counts=driver_counts, SNV_counts=SNV_counts)


class RandomHybridLibrary(BaseHybridLibrary):
    """The hybrid selection model with randomly generated selection coefficients.

    The hybrid model extends the region model to consider SNVs as well. The idea is that SNVs can disrupt the function of genes, thereby altering the selective effect of that gene on the cell's fitness. Here we distinguish between driver SNVs, which have this effect, and passenger SNVs, which have no effect (synonymous). For each region on each chromosome :math:`(k,i)`, randomly assign a second coefficient :math:`\\lambda_{k,i}` in the range (-SNV_coeff, SNV_coeff), where :math:`\\lambda_{k,i} > 0` if :math:`(k,i)` is a TSG but :math:`\\lambda_{k,i} > 0` or :math:`\\lambda_{k,i} < 0` if :math:`(k,i)` is an OG. Additionally, a fraction of NEU regions are set to be essential genes (EG), where if :math:`(k,i)` is an EG then :math:`\\delta_{k,i} = 0` and :math:`\\lambda_{k,i} > 0`, while :math:`\\lambda_{k,i} = 0` for the remaining neutral genes. Fitness is computed as

    .. math::
    
       s_a(c) = \\prod_{k=1}^K \\prod_{i=1}^{m_k} \\Big(1 + \\delta_{k,i}\\cdot(e^{-\\hat{y}_{k,i} \\lambda_{k,i}})\\Big)^{x_{k,i}/p_c} \\cdot \\Big(e^{-\\hat{y}_{k,i} \\lambda_{k,i}}\\Big)^{h(k,i)},

    where :math:`\\hat{y}_{k,i}` is the average number of driver SNVs on any copy of region :math:`i' of chromosome :math:`k' and :math:`h` is an indicator function with :math:`h(k,i) = 1` if region :math:`(k,i)` is an EG and :math:`h(k,i) = 0` otherwise.

    """
    def initialize(
        self, 
        params: Optional[Parameters] = None, 
        CN_coeff: Optional[float] = None,
        SNV_coeff: Optional[float] = None,
        OG_r: Optional[float] = None, 
        TSG_r: Optional[float] = None, 
        EG_r: Optional[float] = None
    ):
        """Method used to initialize selection coefficients. See :ref:`Parameters <parameters>` for an explanation of the parameters.

        Args:
            params (Parameters, optional):
            CN_coeff (float, optional):
            SNV_coeff (float, optional):
            OG_r (float, optional):
            TSG_r (float, optional):
            EG_r (float, optional):
        """
        params = fill_params(params, CN_coeff=CN_coeff, SNV_coeff=SNV_coeff, OG_r=OG_r, TSG_r=TSG_r, EG_r=EG_r)

        NEU_r = 1 - params.OG_r - params.TSG_r - params.EG_r
        for chrname, num_regions in self.regions.items():
            labels = np.random.choice([0, 1, 2, 3], size=num_regions, p=[params.OG_r, params.TSG_r, params.EG_r, NEU_r])
            for r,l in enumerate(labels):
                if l == 0:
                    # OG
                    self.delta[chrname].append(np.random.uniform(0, params.CN_coeff))
                    self.lam[chrname].append(params.SNV_coeff*((-1)**np.random.randint(0, 2)))
                    self.drivers[chrname].append(r)
                elif l == 1:
                    # TSG
                    self.delta[chrname].append(np.random.uniform(-params.CN_coeff, 0))
                    self.lam[chrname].append(params.SNV_coeff)
                    self.drivers[chrname].append(r)
                elif l == 2:
                    # EG
                    self.delta[chrname].append(0)
                    self.lam[chrname].append(params.SNV_coeff)
                    self.essential[chrname].append(r)
                else:
                    # NEU
                    self.lam[chrname].append(0)
                    self.delta[chrname].append(0)
        self.max_distinct_driv = int(self.max_distinct_driv_ratio * sum([len(self.drivers[chrname]) for chrname in self.drivers]))
        self.init_base_fit()
        self.init_max_fit()

class FittedHybridLibrary(BaseHybridLibrary):
    """Uses the same model as :code:`RandomHybridLibrary`, but the selection coefficients are given by the user.

    """
    def initialize(
        self, 
        filepath: Optional[str] = None,
        delta: Optional[Dict] = None,
        lam: Optional[Dict] = None,
    ):
        """Method used to initialize selection coefficients.

        Args:
            filepath (str, optional): Path to a tsv file containing the selection coefficients. The first entry is the chromosome name, followed by the selection coefficient pairs in order of reference regions. Each pair is of the form x,y where x is the CN coefficient and y is the SNV coefficient (:math:`\\delta` and :math:`\\lambda`).
            delta (dict, optional): Instead of passing a file, users can pass a dictionary where keys are chromosome names and values are an :math:`m_k`-length array of CN coefficients. Must also be passed with lam in the same form.
            lam (dict, optional): Instead of passing a file, users can pass a dictionary where keys are chromosome names and values are an :math:`m_k`-length array of SNV coefficients. Must also be passed with delta in the same form.
        
        """
        if delta is not None and lam is not None:
            if set(delta.keys()) != set(lam.keys()):
                raise ValueError("Provided delta and lam inputs should use the same chromosome names.")
            if len(set(delta.keys()) - set(self.chrom_lens.keys())) > 0:
                raise ValueError("Provided selection coefficient dictionary contains unrecognized chromosomes.")
            
            deltas, lam = {}, {}
            for chrname,v in delta.items():
                if not isinstance(v, (list, tuple, np.ndarray)) or type(v) != type(lam[chrname]):
                    raise ValueError("Provided selection coefficient dictionary values (both delta and lam) must be an array or list containing float values.")
                if len(v) != len(lam[chrname]):
                    raise ValueError(f"Number of deltas for {chrname} is different to the number of lambdas.")
                if len(v) > self.regions[chrname]:
                    raise ValueError(f"Number of deltas for {chrname} exceeds the number of regions specified by chrom_lens and region_len on library initialization.")
                for i in range(len(v)):
                    if not isinstance(v[i], (float, np.floating)) or not isinstance(lam[chrname][i], (float, np.floating)):
                        raise ValueError(f"Provided deltas are neither floats nor np.floats.")
                
                v1 = [float(x) for x in v]
                v2 = [float(x) for x in lam[chrname]]
                while len(v1) != self.regions[chrname]:
                    v1.append(0)
                    v2.append(0)
                deltas[chrname] = v1
                lam[chrname] = v2
        
        elif filepath is not None:
            assert os.path.isfile(filepath), 'File does not exist at specified path.'
            deltas, lam = {}, {}
            with open(filepath) as f:
                for line in f:
                    info = line.strip().split('\t')
                    chrname = info[0]
                    if chrname not in self.chrom_lens:
                        raise ValueError("Provided arm selection coefficients file contains unrecognized chromosome.")
                    deltas[chrname] = []
                    lam[chrname] = []
                    for x in info[1:]:
                        y = x.split(',')
                        deltas[chrname].append(float(y[0]))
                        lam[chrname].append(float(y[1]))

        else:
            raise ValueError("Must pass either filepath or delta and lam.")
        
        for chrname in self.regions:
            if chrname not in deltas:
                deltas[chrname] = [0 for _ in range(self.regions[chrname])]
                lam[chrname] = [0 for _ in range(self.regions[chrname])]
        self.delta = deltas
        self.lam = lam
        
        for chrname,v in deltas.items():
            for i,x in enumerate(v):
                if x != 0:
                    self.drivers[chrname].append(i)
                elif x == 0 and lam[chrname][i] != 0:
                    self.essential[chrname].append(i)

        self.max_distinct_driv = int(self.max_distinct_driv_ratio * sum([len(self.drivers[chrname]) for chrname in self.drivers]))
        self.init_base_fit()
        self.init_max_fit()