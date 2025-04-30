import numpy as np

from sistem.anatomy.base_anatomy import BaseAnatomy

class StaticAnatomy(BaseAnatomy):
    def __init__(self, **kargs):
        super().__init__(**kargs)
        self.tmatrix = None
    
    def initialize_distances(self, method='random', matrix=None, path=None):
        super().initialize_distances(method=method, matrix=matrix, path=path)
        
        rates = [(1/x)*self.eps if x != 0 else 0 for x in self.dists[:, 0]]
        idxs = np.triu_indices(self.nsites)
        self.tmatrix = np.zeros((self.nsites, self.nsites))
        self.tmatrix[idxs] = rates
        self.tmatrix = self.tmatrix + self.tmatrix.T - np.diag(np.diag(self.tmatrix))
        self.initialized = True

    def get_oneway_transition_probabilities(self, clone):
        return self.tmatrix[clone.site]
    
    def get_total_transition_probability(self, clone):
        sprobs = self.tmatrix[clone.site]
        tot_prob = 1 - np.prod([1 - x for x in sprobs])
        return tot_prob