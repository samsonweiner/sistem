import numpy as np

from sistem.anatomy.base_anatomy import BaseAnatomy

class GenotypeAnatomy(BaseAnatomy):
    def __init__(self, **kargs):
        super().__init__(**kargs)
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