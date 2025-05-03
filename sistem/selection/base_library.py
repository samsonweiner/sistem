from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional, Union, Dict
from abc import ABC, abstractmethod

from sistem.genome.utils import get_num_regions, fixed_chrom_lens
from sistem.parameters import Parameters, fill_params
from sistem.utilities.utilities import init_mbit, get_reg_id

@dataclass
class Attributes:
    driver_counts: Optional[Dict] = field(default=None)
    arm_counts: Optional[Dict] = field(default=None)
    SNV_counts: Optional[Dict] = field(default=None)
    blacklisted_regions: Optional[Dict] = field(default=None)

class BaseLibrary(ABC):
    """Abstract base class for selection libraries. If creating a custom selection library, users must create a child class which inherits from this class. See :ref:`Parameters <parameters>` for an explanation of the constructor parameters.

    Args:
        params (Parameters, optional):
        chrom_lens (dict, optional):
        arm_ratios (float, dict, optional):
        region_len (int, optional):
        max_distinct_driv_ratio (float, optional):
        max_region_CN (int, optional):
        max_region_SNV (int, optional):
        max_ploidy (int, float, optional):
        min_ploidy (int, float, optional):

    Attributes:
        is_driver_region_model (bool): True if individual regions can be drivers, False otherwise.
        is_driver_SNV_model (bool): True if individual SNVs can impact fitness, False otherwise.

    """
    is_driver_region_model = None
    is_driver_SNV_model = None

    def __init__(
        self, 
        params: Optional[Parameters] = None,
        chrom_lens: Optional[dict] = None, 
        arm_ratios: Optional[Union[float, dict]] = None, 
        region_len: Optional[int] = None, 
        max_distinct_driv_ratio: Optional[float] = None, 
        max_region_CN: Optional[int] = None, 
        max_region_SNV: Optional[int] = None, 
        max_ploidy: Optional[Union[int, float]] = None, 
        min_ploidy: Optional[Union[int, float]] = None, 
    ):
        params = fill_params(params, chrom_lens=chrom_lens, arm_ratios=arm_ratios, region_len=region_len, max_distinct_driv_ratio=max_distinct_driv_ratio, max_region_CN=max_region_CN, max_region_SNV=max_region_SNV, max_ploidy=max_ploidy, min_ploidy=min_ploidy)
        self.chrom_lens = params.chrom_lens
        self.region_len = params.region_len
        self.arm_ratios = params.arm_ratios
        self.regions = get_num_regions(self.chrom_lens, self.region_len)

        self.max_distinct_driv_ratio = params.max_distinct_driv_ratio
        self.max_distinct_driv = 0
        self.max_region_CN = params.max_region_CN
        self.max_region_SNV = params.max_region_SNV
        self.max_ploidy = params.max_ploidy
        self.min_ploidy = params.min_ploidy

        self.delta = {chrname: [] for chrname in self.regions.keys()}
        self.base_fit = 0
        self.max_fit = 0
        
        self._ndiploid_regions = sum(list(self.regions.values()))
        max_region_id = max(self.regions.values()) 
        self._max_region_id = max_region_id
        init_mbit(max_region_id.bit_length())

    def __getitem__(self, chrname):
        return self.delta[chrname]

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        # Enforce that `is_driver_region_model` is defined and valid
        if not hasattr(cls, 'is_driver_region_model') or not isinstance(cls.is_driver_region_model, bool):
            raise ValueError(
                f"Class {cls.__name__} must define a 'is_driver_region_model' attribute set to True or False, depending on whether or not the model involves driver regions."
            )

        # Enforce that `is_driver_SNV_model` is defined and valid
        if not hasattr(cls, 'is_driver_SNV_model') or not isinstance(cls.is_driver_SNV_model, bool):
            raise ValueError(
                f"Class {cls.__name__} must define a 'is_driver_SNV_model' attribute set to True or False, depending on whether or not the model involves driver SNVs."
            )

    @abstractmethod
    def compute_fitness(self, clone):
        """Computes the fitness of a clone.
        
        """

    @abstractmethod
    def init_base_fit(self):
        """Sets the base_fit attribute to the fitness of an unmutated Cell.
        """
        pass

    @abstractmethod
    def init_max_fit(self):
        """Sets the max_fit attribute to an estimate of the maximum possible fitness.
        
        """
        pass

    @abstractmethod
    def check_viability(self, clone):
        """Checks if clone passes viability checkpoints based on mutated driver stats. 

        Args:
            clone (Clone): The clone object in question.

        Returns:
            (bool): True if passes, False if not.
        
        """
        pass

    @abstractmethod
    def update_stats(self, clone, chromosome, start, end, mag=1):
        """Function which, upon a mutation occurring, updates intermediate stats stores in the clones Attributes object to make fitness computation more efficient.

        Arguments:
            clone (Clone): The clone undergoing a mutation.
            chromosome (Chromosome): The Chromosome object being mutated.
            start (int): The index of the starting region.
            end (int): The index of the ending region.
            mag (int): If the event is a CNA (start != end), then it is the number of copies being gained. If the event is an SNV (start == end), then it interprets an SNV being added in region start.
        """
        pass

    @abstractmethod
    def init_attributes(self, clone):
        """Creates an Attributes object for a clone and updates it to reflect its genome state.

        Args:
            clone (Clone): The clone object in question.
        """

    @abstractmethod
    def get_driver_start_regions(self, cell, chromosome, size):
        """Identifies the possible start region indices of a driver event.

        Args:
            cell (Cell): The current cell.
            chromosome (Chromosome): The current chromosome.
            size (int): The size in number of regions of the event (1 for SNVs, >=1 for CNAs).
        """

    @abstractmethod
    def get_passenger_start_regions(self, cell, chromosome, size):
        """Identifies the possible start region indices of a passenger event.

        Args:
            cell (Cell): The current cell.
            chromosome (Chromosome): The current chromosome.
            size (int): The size in number of regions of the event (1 for SNVs, >=1 for CNAs).
        """

    @abstractmethod
    def initialize(self, **kwargs):
        """Function used to initialize the selection coefficients.
        
        """
        pass