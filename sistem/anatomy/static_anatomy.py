import numpy as np
from typing import Optional, Union, List

from sistem.anatomy.base_anatomy import BaseAnatomy
from sistem.selection import BaseLibrary

class StaticAnatomy(BaseAnatomy):
    """The static migration model.
    
    In this model, the probabiltiy of migrations from site :math:`a` to site :math:`b` depend on the distance :math:`d(a,b)`, but are the same for each cell in the site at each generation. 
    """
    def __init__(self, libraries: Optional[Union[BaseLibrary, List[BaseLibrary]]], **kargs):
        super().__init__(libraries=libraries, **kargs)
        self.tmatrix = None
    
    def initialize_distances(self, method='random', matrix=None, path=None):
        super().initialize_distances(method=method, matrix=matrix, path=path)
        
        rates = [(1/x)*self.eps if x != 0 else 0 for x in self.dists[:, 0]]
        idxs = np.triu_indices(self.nsites)
        self._tmatrix = np.zeros((self.nsites, self.nsites))
        self._tmatrix[idxs] = rates
        self._tmatrix = self._tmatrix + self._tmatrix.T - np.diag(np.diag(self._tmatrix))
        self.initialized = True

    def get_oneway_transition_probabilities(self, clone):
        return self._tmatrix[clone.site]
    
    def get_total_transition_probability(self, clone):
        sprobs = self._tmatrix[clone.site]
        tot_prob = 1 - np.prod([1 - x for x in sprobs])
        return tot_prob