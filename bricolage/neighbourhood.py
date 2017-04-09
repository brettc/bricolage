from bricolage.core_ext import Collection
from math import log as mathln
import numpy as np

class NetworkNeighbourhood(object):
    """Gather a sample of nearby mutants

    Samples can be simple a one mutation away (exactly what this means depends
    on the model of mutation), or they can be a number of mutations away,
    where the number of steps is drawn from poisson distribution, ensuring
    that the expected number of single-step walks is proportion of one-step
    neighbourhood.
    """
    def __init__(self, center, sample_size,
                 one_step_proportion=1.0,
                 lamb=0.0,
                 steps=1,
                 ):
        """There three ways to generation a neighbourhood.
            - specify the number of steps
            - specifiy the proportion that will be one step.
            - specify a lambda parameter / mean for sampling usig poisson.

        """
        self.sample_size = sample_size
        self.center = center

        if one_step_proportion < 1.0 or lamb > 0.0:
            if lamb > 0.0:
                self.lamb = lamb
            else:
                # Set lambda so we get the right proportion of one-step neighbours
                self.lamb = -mathln(one_step_proportion)

            # We want at least one mutation, so add one.
            self.mutations = np.random.poisson(self.lamb, self.sample_size) + 1
        else:
            self.mutations = np.ones(self.sample_size, dtype=int) * steps
            # print self.mutations[:10]


        self.neighbours = Collection(center.factory)
        self.neighbours.fill_with_mutations(center, self.mutations)


class PopulationNeighbourhood(object):
    """The same thing, for populations"""
    def __init__(self, population, sample_size,
                 one_step_proportion=1.0,
                 lamb=0.0,
                 steps=1,
                 ):
        """See above"""

        self.one_step_proportion = one_step_proportion
        self.sample_size = sample_size
        self.population = population
        self.neighbours = Collection(population.factory)

        if one_step_proportion < 1.0 or lamb > 0.0:
            if lamb > 0.0:
                self.lamb = lamb
            else:
                # Set lambda so we get the right proportion of one-step neighbours
                self.lamb = -mathln(one_step_proportion)

        # We simple iterate over all networks and add to our neighbourhood
        for n in population:
            if lamb > 0.0:
                mutes = np.random.poisson(self.lamb, self.sample_size) + 1
            else:
                mutes = np.ones(self.sample_size, dtype=int) * steps

            network_neighbours = Collection(population.factory)
            network_neighbours.fill_with_mutations(n, mutes)

            # We just add to the main population
            self.neighbours.extend(network_neighbours)
