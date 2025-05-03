import numpy as np
from typing import Optional, Union, List

from sistem.anatomy.base_anatomy import BaseAnatomy
from sistem.selection import BaseLibrary

class GenotypeAnatomy(BaseAnatomy):
    """The genotype migration model.
    
    In this model, the probability of cell :math:`c` migrating from site :math:`a` to site :math:`b` increases as :math:`c` gains beneficial mutations with respect to the selection landscape of :math:`b`. In particular, distance is defined as

    .. math::

        d(a,b) = \\Big(\\log{\\frac{s_b(c)}{\\hat{s}_{b}}}\\Big)^{-1},

    where :math:`\\hat{s}_{b}` is the fitness of a non-mutated cell in site :math:`b`. In other words, the probability that cell :math:`c` successfully migrates to :math:`a` site :math:`b` increases as the cell acquires beneficial mutations. When used in conjunction with site-specific selection libraries, migration probabilities can also reflect compatibility with the the fitness landscape :math:`b`.

    """
    def __init__(self, libraries: Optional[Union[BaseLibrary, List[BaseLibrary]]], **kargs):
        super().__init__(libraries=libraries, **kargs)
        self.initialized = True

    def get_oneway_transition_probabilities(self, clone):
        s = clone.site
        probs = []
        for i in range(len(self.libraries)):
            if i == s:
                probs.append(0)
            else:
                multiplier = max(np.log(self.libraries[i].compute_fitness(clone) / self.libraries[i].base_fit), 1)
                probs.append(multiplier * self.eps)
        return np.array(probs)
    
    def get_total_transition_probability(self, clone):
        sprobs = self.get_oneway_transition_probabilities(clone)
        tot_prob = 1 - np.prod([1 - x for x in sprobs])
        return tot_prob