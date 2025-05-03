import numpy as np

import string
from itertools import product
from abc import ABC, abstractmethod
from typing import Optional, Union, List

from sistem.selection import BaseLibrary
from sistem.anatomy.utils import create_random_distances
from sistem.utilities.IO import load_distances_from_file
from sistem.parameters import Parameters, fill_params

class BaseAnatomy(ABC):
    """Abstract base class for migration models. If creating a custom migration model, users must create a child class which inherits from this class. See :ref:`Parameters <parameters>` for an explanation of the constructor parameters.

    Migration are used to define a distance function :math:`d(a,b)` between each pair of sites a and b. The probability that cell :math:`c` migrates from :math:`a` to :math:`b` is given as 

    .. math::

        \\alpha_c(a, b) = \\left\\{\\begin{array}{lr} 0 & \\text{if } a = b, \\\\ \\frac{1}{d(a,b)} \\cdot \\epsilon & \\text{otherwise}\\end{array}\\right\\},

    where :math:`\\epsilon` is a per-generation baseline migration rate.

    Args:
        libraries (BaseLibrary, List[BaseLibrary]): The primary site selection library, or a list of libraries (one for each site). Required.
        params (Parameters, optional):
        nsites (int, optional):
        growth_rate (float, optional):
        max_growth_rate_multiplier (int, float, optional):
        capacities (int, List[int], optional):
        epsilon (float, optional):
        t_max (int, optional):
        N0 (int, optional):

    Attributes:
        exp_pop (List[list]): A list of lists, one for each site, containing the expected population size for each generation t=0,...,t_max.
        dists (array): The pairwise distances between anatomical sites, in order of the upper triangular indices.
        points (array): If initialized with random distances, contains the 2d coordinates of each site used to derive the distances.
    """
    def __init__(
        self, 
        libraries: Optional[Union[BaseLibrary, List[BaseLibrary]]] = None,
        params: Optional[Parameters] = None,
        *,
        nsites: Optional[int] = None, 
        growth_rate: Optional[float] = None, 
        max_growth_rate_multiplier: Optional[Union[int, float]] = None, 
        capacities: Optional[Union[int, List[int]]] = None, 
        epsilon: Optional[float] = None, 
        t_max: Optional[int] = None, 
        N0: Optional[int] = None
    ):
        params = fill_params(params, nsites=nsites, growth_rate=growth_rate, max_growth_rate_multiplier=max_growth_rate_multiplier, capacities=capacities, epsilon=epsilon, t_max=t_max, N0=N0)
        
        self.nsites = params.nsites
        self.growth_rate = params.growth_rate
        self.eps = params.epsilon
        self.t_max = params.t_max
        self.N0 = params.N0
        if params.max_growth_rate_multiplier is None or params.max_growth_rate_multiplier < 1:
            self.max_meta_growth_rate = self.growth_rate
        else:
            self.max_meta_growth_rate = self.growth_rate*params.max_growth_rate_multiplier
        self.capacities = self._check_capacities(params.capacities)
        self.exp_pop = [[] for i in range(self.nsites)]

        chars = list(string.ascii_uppercase)
        chars.remove('P')
        self.site_ids = ['P']
        if self.nsites > 26:
            chars = [''.join(x) for x in product(chars, repeat=2)]
        self.site_ids += chars[:self.nsites - 1]

        self.init_pop_dynamic(0, 0, params.t_max, N0=params.N0)

        self.libraries = self._check_libraries(libraries)
        self.dists = None
        self.points = None
        self.initialized = False

    def _check_capacities(self, capacities):
        if isinstance(capacities, list):
            if len(capacities) == 1:
                return capacities*self.nsites
            if len(capacities) != self.nsites:
                raise ValueError("Number of capacities does not match number of sites.")
            return capacities
        elif isinstance(capacities, (int, float)):
            return [int(capacities)]*self.nsites
        else:
            raise ValueError("The provided 'capacities' argument must be an int or list of ints.")
    
    def _check_libraries(self, libraries):
        if libraries is None:
            raise ValueError("Argument 'libraries' is required and cannot be None. Provide either a single Library class instance or a list of Library class instances.")
        
        if isinstance(libraries, list):
            if not all(isinstance(lib, BaseLibrary) for lib in libraries):
                raise TypeError("The provided 'libraries' list contains elements which are not Library class instances.")
            
            if len(libraries) != self.nsites:
                raise ValueError("If passing a list of libraries, ensure the number of libraries is equal to the number of sites.")
            
            primary_library_type = type(libraries[0])
            if not all(type(lib) is primary_library_type for lib in libraries):
                raise TypeError("All elements must be instances of the same Library class.")
            return libraries

        elif isinstance(libraries, BaseLibrary):
            return [libraries] * self.nsites
        
        else:
            raise TypeError("The provided 'libraries' argument must be a Library class instance or a list of Library class instances.")
    
    def init_pop_dynamic(self, site, t, t_max, N0=None, cur_fit=None, library=None):
        pops = [0 for i in range(t)]
        if N0 is not None:
            pops.append(N0)
            gr = self.growth_rate
        else:
            pops.append(self.N0)
            if self.growth_rate == self.max_meta_growth_rate or library == None:
                gr = self.growth_rate
            else:
                r = np.log(1 + cur_fit - library.base_fit) / np.log(1 + library.max_fit - library.base_fit)
                gr = self.growth_rate + r * (self.max_meta_growth_rate - self.growth_rate)
        for i in range(t, t_max+1):
            pn = pops[i]
            pn1 = pn + gr*pn*(1 - pn/self.capacities[site])
            pops.append(pn1)
        self.exp_pop[site] = pops

    def initialize_distances(self, method='random', matrix=None, path=None):
        """Method used to initialize the organotropism priors. 

        Args:
            method (str): The method used to initialize the distances. Either 'random' or 'precomputed'. If precomputed, must provide either a matrix or path. Note that if precomputed, the distance between a and itself must be 0.
            matrix (list, np.ndarray): A pairwise distance matrix. Must either be square or 1d containing just the upper triangular indices in order. 
            path (str): Instead of passing a python object, users can also provide a path to a file containing the distance matrix. Each row of the matrix should be on a separate line, and values should be separated by 0s. 
        """
        if method == 'random':
            points, norm_dists = create_random_distances(self.nsites)
            self.points = points
            self.dists = norm_dists
        elif method == 'precomputed':
            if matrix is not None:
                if not isinstance(matrix, (list, np.ndarray)):
                    raise ValueError("If providing a precomputed distance matrix, must be either a list or numpy array.")
                matrix = np.array(matrix)
                if matrix.ndim == 1:
                    self.dists = matrix
                else:
                    idxs = np.triu_indices(self.nsites)
                    self.dists = matrix[idxs]
            elif path is not None:
                self.dists = load_distances_from_file(self.nsites, path)
            else:
                raise ValueError("If using precomputed distances, either 'matrix' (list or numpy array) or 'path' (path to file containing distances) must be specified.")

        else:
            raise ValueError("Argument 'method' must be either 'random' or 'precomputed'.")
    
    def create_random_metastatic_libraries(
        self, 
        method: str = 'random',
        params: Optional[Parameters] = None,
        alter_prop: Optional[float] = None, 
        CN_coeff: Optional[float] = None
    ):
        """Creates site-specific selection libraries. The anatomy class defaults to using the same selection library for all sites, and this method is only used if site-specific selection is desired. There are two options for the *method* parameter: 'random' or 'distance'. If random, the number of selection coefficients in each metastatic site will differ from the primary site with a ratio near *alter_prop*. If distance, the number of selection coefficients that differ between all pairs of sites will roughly align with the pairwise distance matrix created with the :code:`initialize_distance` method See :ref:`Parameters <parameters>` for an explanation of the parameters.

        Args:
            method (str): The method used to create the metastatic libraries. Must be either 'random' or 'distance'.
            params (Parameters, optional):
            alter_prop (float, optional):
            CN_coeff (float, optional):
            
        """
        params = fill_params(params, alter_prop=alter_prop, CN_coeff=CN_coeff)

        if self.nsites == 1:
            raise ValueError("Number of sites must be >1 to create metastatic libraries.")
        primary_library = self.libraries[0]
        primary_library_type = type(primary_library)
        if hasattr(primary_library_type, 'create_metastatic_libraries') and callable(getattr(primary_library_type, 'create_metastatic_libraries')):
            self.libraries = primary_library.create_metastatic_libraries(self.nsites, params.alter_prop, params.CN_coeff, method=method, dists=self.dists)
        else:
            raise AttributeError(f"{primary_library_type.__name__} must have a callable 'create_metastatic_libraries' method in order to call this function.")

    
    @abstractmethod
    def get_oneway_transition_probabilities(self, clone):
        """Method which returns an array of the migation probabilities between the clone and all sites at the current generation

        Args:
            clone (Clone, optional): The clone.

        Returns:
            list: The migration probabilities.
        """
        pass

    @abstractmethod
    def get_total_transition_probability(self, clone):
        """Method which returns the total probability that a clone migrates to any other site.

        Args:
            clone (Clone, optional): The clone.

        Returns:
            float: The total migration probability.
        """
        pass

class SimpleAnatomy(BaseAnatomy):
    """The simplest migration model. Use for simulating single-site experiments, or if you want an identical fixed migration probability to each site for every cell at every generation (:math:`d(a,b)=1` for all :math:`a,b`).

    """
    def __init__(self, libraries: Optional[Union[BaseLibrary, List[BaseLibrary]]], **kwargs):
        super().__init__(libraries=libraries, **kwargs)
        self._tmatrix = np.array([[self.eps for j in range(i)] + [0] + [self.eps for j in range(i+1, self.nsites)] for i in range(self.nsites)])
        self.initialized = True

    def get_oneway_transition_probabilities(self, clone):
        return self._tmatrix[clone.site]

    def get_total_transition_probability(self, clone):
        return 1 - (1 - self.eps)**self.nsites