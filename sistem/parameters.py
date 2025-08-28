import os
import warnings
from dataclasses import dataclass, field
from typing import Optional, Union, List, Dict

from sistem.genome import hg38_chrom_lengths_from_cytoband, get_chrom_lens_from_reference

@dataclass
class Parameters:
    """
    Global parameter values for a simulation experiment.

    Attributes:
        N0 (int): The initial number of cells at time t=0 and upon site seeding. **Default: 10.**
        growth_rate (float): The primary site exponential growth rate. **Default: 0.0051.**
        max_growth_rate_multiplier (int, float): This parameter is used to increase the exponential growth rates in metastatic sites. Metastatic growth rates are at minimum *growth_rate* but are increased based on the fitness of the cell which initiates seeding, up to a maximum of *growth_rate* * *max_growth_rate_multiplier*. **Default: 5.**
        capacities (int, list): The carrying capacity (max number of cells) of each anatomical site. If an integer is provided, all sites have that same carrying capacity. Also accepts a list of carrying capacities, one for each site. **Default: 1e7.**
        epsilon (float): A per-cell per-generation baseline migration probability to a site. **Default: 1e-8.**
        lifespan_mean (int, float): The mean cell lifespan in number of generations drawn from an exponential distribution. Corresponds to 1/r, where r is the rate parameter of the exponential. A value of 1 will fix the lifespan to 1 generation for all cells. 

        ref (str, optional): Path to an input reference genome in fasta format. Can be used to initialize chromosome sizes, but is required for generating synthetic sequencing reads. **Default: None.**
        alt_ref (str, optional): Path to an optionally alternate reference genome in fasta format. Allele 0 utilizes ref, while Allele 1 utilizes alt_ref. Use if the goal is to generate allele-specific synthetic sequencing reads. **Default: None.**
        chrom_names (list, optional): A list of chromosome names to use in the reference sequence(s). Only necessary to specify if the reference genome contains superfluous chromosomes which should not be included in the genome. **Default: None.**
        chrom_lens (dict, int, optional): A dictionary describing the size of the genome where keys are chromosome names and values are chromosome lengths in number of base pairs. If kept as None, will automatically populate with chr1-chr22 and lengths derived from the hg38 human reference genome. If ref is specified, will automatically populate based on the provided sequences. Additionally, if an integer is passed, then all chromosomes 1-*num_chrom* (see below) will have the same given length. **Default: None.**
        num_chrom (int): The number of chromosomes in the genome. Use only if an int is passed to *chrom_lens*, or if you want to use the first *num_chrom* human reference chromosomes if *chrom_lens* is set to None. Otherwise, will update accordingly. **Default: 22.**
        region_len (int): SISTEM utilizes a simplified genome representation whereby chromosome sequences are partitioned into non-overlapping regions of uniform size *region_len* (in base pairs). Higher values reduce memory burden, while lower values increase simulation resolution. **Default: 5e6.**

        arm_ratios (float, dict): The ratio of the small chromosome arm length to the total chrosome length. If chrom_lens is None, will utilize arm ratios derived from the hg38 human reference genome. Can pass a dictionary where keys are chromosome names and values are ratios, or a single ratio used by all chromosomes. **Default: 0.5.**

        focal_driver_rate (float): The probability of acquiring a driver focal (segmental) CNA at each generation. **Default: 1e-4.**
        focal_pass_rate (float): The probability of acquiring a passenger focal (segmental) CNA at each generation. **Default: = 0.01.**
        length_mean (int, float): The mean number of regions a focal CNA spans. Length is drawn from an exponential distribution. **Default = 1.5.**
        focal_gain_rate (float): The probability that a focal CNA is amplification (gain) versus a deletion (loss) **Default = 0.5.**
        mag_mean (int, float): The mean number of additional copies gained during a focal amplification CNA. Amplification magnitude is drawn from a geometric distribution. Default = 1.2.**
        SNV_driver_rate (float): The probability of acquiring a driver SNV at each generation. **Default: 1e-4**
        SNV_pass_rate (float): The probability of acquiring a passenger SNV at each generation. **Default: 0.01.**
        arm_rate (float): The probability of acquiring a chromosome-arm CNA at each generation. **Default: 1e-5.**
        chromosomal_rate (float): The probability of acquiring a whole-chromosomal CNA at each generation. **Default: 1e-6.**
        chrom_dup_rate (float): The probability that a chromosome-arm CNA or a whole-chromosomal CNA is a duplication versus a deletion. **Default: 1e-5.**
        WGD_rate (float): The probability of acquiring a WGD at each generation. **Default: 1e-8.**

        CN_coeff (float): The maximum CN selection coefficient magnitude. Used only for random initialization. **Default: 0.25.**
        SNV_coeff (float): The maximum region SNV selection coefficient magnitude. Used only for random initialization. **Default: 0.1.**
        OG_r (float): The ratio of regions which are OGs. Used only for random initialization in the Region/Hybrid Selection Model. **Default: 0.05.**
        TSG_r (float): The ratio of regions which are TSGs. Used only for random initialization in the Region/Hybrid Selection Model.. **Default: 0.05.**
        EG_r (float): The ratio of regions which are EGs (essential genes). Used only for random initialization in the Hybrid Selection Model. **Default: 0.05.**
        alter_prop (float): Parameter for creating site-specific metastatic libraries. Represents the fraction of driver selection coefficients to alter in each site if method is 'random', or of the farthest site if method is 'distance', with the number of altered coefficients scaled accordingly for the rest. **Default: 0.3.**

        max_region_CN (int): Viability checkpoint parameter. Represents the maximum CN of a driver region. **Default: 10.**
        max_region_SNV (int): Viability checkpoint parameter. Represents the maximum number of driver SNVs in a single region. **Default: 10.**
        max_ploidy (int, float): Viability checkpoint parameter. Represents the maximum ploidy of a cell. **Default: 8.**
        min_ploidy (int, float): Viability checkpoint parameter. Represents the minimum ploidy of a cell. **Default: 1.5.**
        max_distinct_driv_ratio (float): Viability checkpoint parameter. Represents the maximum number of distinct drivers which can be mutated. **Default: 0.8.**

        t_max (int): The maximum number of generations to run. **Default: 6000.**
        min_detectable (int): Terminates the growth simulator when the number of cells present in each anatomical site is atleast *min_detectable*. **Default: 5e-5.**
        ncells_prim (int): The number of cells to sample from the primary site. **Default: 100.**
        ncells_meta (int): The number of cells to sample from the metastatic sites. **Default: 100.**
        ncells_normal (int): The number of normal cells to diluate the primary site with. If generating clonal lineages, a relative number of normal cells will be sampled from the metastatic sites as well. Passing 1 will add a convenient normal cell outgroup to the simulated lineage tree. **Default: 1.**
        min_mut_fraction (float): In SISTEM, clones are defined by cells with a unique sequence of driver mutations, but this loose definition means that distinct 'clones' appearing in the clonal lineage tree may differ by a just a few small mutations. The *min_mut_fraction* parameter can help make the clones present in the tree more distinct. It describes the minimum frequency a clone's genotype must occur in the sampled cells to remain in the tree. If possible, multiple clones in a subtree will merge together with a common genotype to remain above *min_mut_fraction*, otherwise they are pruned. **Default: 0.05**.
 
        bin_size (int, optional): The size in base pairs of the copy number windows/segments in the final output profiles. Essentially groups consecutive regions together into larger bins and computes the mean. If kept as None, the bin size will be set equal to the *region_len*. This will be desirable for the majority of cases. Only specify if simulating with a region length that is smaller than practical (e.g. if *region_len* is less than 100kbp for single-cell data). **Default: None.**
        coverage (int, float): The average number of reads which cover any given base pair in the genome. When used to generate single-cell read counts or DNA-seq reads, it is recommended to use a low value (<0.2), whereas if used to generate bulk read counts, it is recommended to use a high value (>50). **Default: 0.1.**
        read_len (int): The length of the paired-end reads. Together with coverage, used to compute expected read counts. **Default: 150.**
        seq_error (float): Per-base pair sequencing error rate. Only used for generating scDNA-seq reads **Default: 0.02.**
        lorenz_y (float): Used to introduce coverage non-uniformity when generating single-cell read counts and DNA-seq. In a nutshell, coverage uniformity is parameterized by a point on the lorenz curve (0.5, *lorenz_y*). Default value of 0.5 means maximally uniform, and decreasing *lorenz_y* down to 0 will decrease uniformity. Only specify if evaluating conditions under non-uniform coverage. Be aware that this parameter is extremely sensitive, and values <=0.4 will lead to highly non-uniform distributions. **Default: 0.5.**

        out_dir (str): The path to the output directly. **Default: 'out'.**
        log_path (str, optional): The path to the log file. By default, will write to *out_dir*/sim.log. **Default: None.**
        num_processors (int): Number of processors to use when generating synthetic scDNA-seq reads. Will not speed up the other steps of the simulator. **Default: 1.**

    """
    N0: int = 10
    growth_rate: float = 0.0051
    max_growth_rate_multiplier: Union[int, float] = 5
    capacities: Union[int, List[int]] = field(default_factory=lambda: [1e7])
    nsites: int = 1
    epsilon: float = 1e-8
    lifespan_mean: Union[int, float] = 1
    
    ref: Optional[str] = None
    alt_ref: Optional[str] = None
    chrom_names: Optional[List[str]] = None
    chrom_lens: Optional[Union[int, Dict]] = None
    num_chroms: int = 22
    region_len: int = field(default=int(5e6))
    arm_ratios: Union[float, Dict] = field(default=0.5)

    focal_driver_rate: float = 5e-4
    SNV_driver_rate: float = 5e-4
    arm_rate: float = 1e-4
    chromosomal_rate: float = 5e-5
    WGD_rate: float = 1e-8
    focal_gain_rate: float = 0.5
    chrom_dup_rate: float = 0.5
    length_mean: Union[int, float] = 1.5
    mag_mean: Union[int, float] = 1.2

    focal_pass_rate: float = 0.01
    SNV_pass_rate: float = 0.01

    # region library
    CN_coeff: float = 0.25
    SNV_coeff: float = 0.1
    OG_r: float = 0.05
    TSG_r: float = 0.05
    EG_r: float = 0.05
    alter_prop: float = 0.3

    max_region_CN: int = 10
    max_region_SNV: int = 10
    max_ploidy: Union[int, float] = 8
    min_ploidy: Union[int, float] = 1.5
    max_distinct_driv_ratio: float = 0.8

    t_max: int = 6000
    min_detectable: int = 5e5
    ncells_prim: int = 100
    ncells_meta: int = 100
    ncells_normal: int = 1
    min_mut_fraction: float = 0.05

    bin_size: Optional[int] = None
    coverage: Union[int, float] = 0.1
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

        self.region_len = int(self.region_len)
        self._process_genome()

        if self.lifespan_mean < 0:
            raise ValueError("Lifespan must be >= 1.")

        if self.CN_coeff < 0 or self.CN_coeff > 1:
            raise ValueError("CN_coeff must be between 0-1.")
        
        if self.SNV_coeff < 0 or self.SNV_coeff > 1:
            raise ValueError("SNV_coeff must be between 0-1.")

        if self.alter_prop < 0 or self.alter_prop > 1:
            raise ValueError("alter_prop must be between 0-1.")

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