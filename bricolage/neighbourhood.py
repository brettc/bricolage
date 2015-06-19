from bricolage.core_ext import Collection
from math import log as mathln
import numpy as np

class NeighbourhoodSample(object):
    """Gather a sample of nearby mutants

    Samples can be simple a one mutation away (exactly what this means depends
    on the model of mutation), or they can be a number of mutations away,
    where the number of steps is drawn from poisson distribution, ensuring
    that the expected number of single-step walks is proportion of one-step
    neighbourhood.
    """
    def __init__(self, center, sample_size, one_step_proportion=1.0):
        # Sample from the neighbourhood
        #
        # @center: our focal network.
        # @one_step_proportion: the proportion of the sample that we want to
        # have only one step mutations
        #
        self.one_step_proportion = one_step_proportion
        self.sample_size = sample_size
        self.center = center

        if one_step_proportion < 1.0:
            # Set lambda so we get the right proportion of one-step neighbours
            lmbda = -mathln(self.one_step_proportion)

            # We want at least one mutation, so add one. Note that the proportion
            # of these that equals 1 should be one_step_proportion.
            self.mutations = np.random.poisson(lmbda, self.sample_size) + 1
        else:
            self.mutations = np.ones(self.sample_size, dtype=int)

        self.neighbours = Collection(center.constructor)
        self.neighbours.fill_with_mutations(center, self.mutations)
        


