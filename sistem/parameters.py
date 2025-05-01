import os
import warnings
from dataclasses import dataclass, field
from typing import Optional, Union, List, Dict

from sistem.genome.utils import hg38_chrom_lengths_from_cytoband, get_chrom_lens_from_reference

@dataclass
class Parameters:
    """
    Global parameter values for a simulation experiment.

    Attributes:
        N0 (int): The initial number of cells at time t=0.
        growth_rate (float): The baseline exponential growth rate.
        max_growth_rate_multiplier: This parameter multiplied by growth_rate represents the maximum exponential growth rate of metastatic tumors.
        capacities (Union[int, List[str]]): The carrying capacity (max number of cells) of each anatomical site. 
    """
    N0: int = 10
    growth_rate: float = 0.005046761847658182
    max_growth_rate_multiplier: Union[int, float] = 5
    capacities: Union[int, List[int]] = field(default_factory=lambda: [1e7])
    nsites: int = 1
    epsilon: float = 1e-8
    
    ref: Optional[str] = None
    alt_ref: Optional[str] = None
    chrom_names: Optional[List[str]] = None
    chrom_lens: Optional[Union[int, Dict]] = None
    num_chroms: int = 22
    region_len: int = field(default=int(1e6))
    arm_ratios: Union[float, Dict] = field(default=0.5)

    focal_driver_rate: float = 1e-4
    SNV_driver_rate: float = 0
    arm_rate: float = 1e-5
    chromosomal_rate: float = 1e-6
    WGD_rate: float = 1e-8
    focal_gain_rate: float = 0.5
    chrom_dup_rate: float = 0.5
    length_mean: Union[int, float] = 1.5
    mag_mean: Union[int, float] = 1.2

    CNA_pass_rate: float = 0.01
    SNV_pass_rate: float = 0.01

    # region library
    CN_coeff: float = 0.1
    SNV_coeff: float = 0.1
    OG_r: float = 0.05
    TSG_r: float = 0.05
    EG_r: float = 0.05
    alter_prop: float = 0.1

    max_region_CN: float = 10
    max_region_SNV: float = 10
    max_ploidy: Union[int, float] = 8
    min_ploidy: Union[int, float] = 1.5
    max_distinct_driv_ratio: float = 0.75

    t_max: int = 6000
    min_detectable: int = 5e5
    ncells_prim: int = 100
    ncells_meta: int = 100
    ncells_normal: int = 1
    min_mut_fraction: float = 0

    bin_size: Optional[int] = None
    coverage: Union[int, float] = 1
    read_len: int = 150
    seq_error: float = 0.02
    lorenz_y: float = 0.5

    out_dir: str = 'out'
    log_path: Optional[str] = None
    num_processors: int = 1

    def __post_init__(self):
        if self.nsites < 1 or self.nsites > 100:
            raise ValueError("Number of sites must be greater than 0 and less than 100.")
        if self.nsites >= 15:
            warnings.warn("Extremely large number of sites. Consider choosing a smaller value.")

        self._process_genome()

        if self.growth_rate <= 0:
            raise ValueError("growth_rate must be greater than 0.")
        if self.max_growth_rate_multiplier < 1:
            raise ValueError("max_growth_rate_multiplier must be 1 or greater.")
        if self.bin_size is None or self.bin_size < self.region_len:
            self.bin_size = self.region_len
        else:
            self.bin_size = int(self.bin_size)
        
        if self.ncells_prim > self.min_detectable or self.ncells_meta > self.min_detectable:
            raise ValueError("min_detectable must be greater than ncells_prim and ncells_meta")
        
        if self.out_dir is None:
            raise ValueError("out_dir must be provided.")

        if self.out_dir is not None:
            if not os.path.isdir(self.out_dir):
                os.makedirs(self.out_dir)
        
        if self.log_path is None:
            self.log_path = os.path.join(self.out_dir, 'sim.log')
    
    def _process_genome(self):
        if self.chrom_names is not None:
            for chrname in self.chrom_names:
                if 'X' in chrname or 'Y' in chrname:
                    warnings.warn("Sex chromosomes detected in provided chrom names, but SISTEM treats all chromosomes as autosomes.")
                    break
            self.num_chroms = len(self.chrom_names)
        if self.num_chroms > 22:
            warnings.warn("Large number of chromosomes detected (>22).")

        if self.ref is not None:
            if not os.path.exists(self.ref):
                raise ValueError("Cannot find file at the provided ref path.")
            chrom_lens = get_chrom_lens_from_reference(self.ref, chrom_names=self.chrom_names)
            if len(chrom_lens) > 22:
                warnings.warn("Provided reference genome contains a large number of chromosomes (>22). Consider specifying which chromosomes to use from the reference genome with the chrom_names parameter. Example: chrom_names = ['chr1', 'chr2', ..., 'chr22'].")
            for chrname in chrom_lens:
                if 'X' in chrname or 'Y' in chrname:
                    warnings.warn("Provided reference genome may contain sex chromosomes (SISTEM treats all chromosomes as autosomes). Consider specifying which chromosomes to use from the reference genome with the chrom_names parameter. Example: chrom_names = ['chr1', 'chr2', ..., 'chr22'].")
            self.chrom_lens = chrom_lens
            if len(chrom_lens) == 22:
                cl, ar = hg38_chrom_lengths_from_cytoband()
                self.arm_ratios = ar

            if self.alt_ref is not None:
                if not os.path.exists(self.alt_ref):
                    raise ValueError("Cannot find file at the provided alt_ref path.")
                alt_chrom_lens = get_chrom_lens_from_reference(self.alt_ref, chrom_names=self.chrom_names)
                if set(chrom_lens.keys()) != set(alt_chrom_lens).keys():
                    raise ValueError("Chromosome names in ref and alt_ref do not match.")
        
        elif self.chrom_lens is None:
            cl, ar = hg38_chrom_lengths_from_cytoband()
            if self.num_chroms > 22:
                raise ValueError("When no reference is passed and chrom_lens is None, hg38 chrom lengths are used. Cannot set num_chroms > 22.")
            if self.num_chroms < 22:
                chrnames = [f'chr{i}' for i in range(1, self.num_chroms + 1)]
                cl = {i: v for i,v in cl.items() if i in chrnames}
                ar = {i: v for i,v in ar.items() if i in chrnames}
            self.chrom_lens = cl
            self.arm_ratios = ar
        
        elif isinstance(self.chrom_lens, int):
            if self.chrom_names is not None:
                self.chrom_lens = {c: self.chrom_lens for c in self.chrom_names}
            else:
                self.chrom_lens = {f'chr{i}': self.chrom_lens for i in range(1, self.num_chroms + 1)}

        elif isinstance(self.chrom_lens, dict):
            self.num_chroms = len(self.chrom_lens)


def fill_params(params, **kwargs):
    """
    Create a Parameters by updating the `params` object based on provided keyword arguments, or creating a new one with the keywords.

    Arguments:
        params: The Parameters object to update.
        **kwargs: Keyword arguments where the keys correspond to the parameter names 
                  to be updated and the values are the new values for those parameters.

    Returns:
        A new Parameters object with updated values.
    """
    # Create a copy of the original params to avoid modifying the original
    if params is None:
        keywords = {}
        for key, value in kwargs.items():
            if value is not None:
                keywords[key] = value
        updated_params = Parameters(**keywords)
        return updated_params
    else:
        params_dict = params.__dict__.copy()
        for key, value in kwargs.items():
            if key in params_dict and value is not None:
                params_dict[key] = value
        updated_params = Parameters(**params_dict)
        return updated_params