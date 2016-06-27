# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API

import cython
import numpy

# cimports
cimport numpy as np
import numpy
from cython.operator import dereference as deref, preincrement as preinc

cdef class BaseTarget:
    def __cinit__(self):
        self._base = NULL

    def __dealloc__(self):
        if self._base != NULL:
            del self._base

    def assess(self, Network net):
        # assert net.factory.world is self.world
        return self._base.assess(deref(net._this));

    def assess_collection(self, CollectionBase coll):
        cdef cRates fitnesses
        self._base.assess_networks(
            deref(coll._collection), fitnesses)

        # Numpy conversion
        # TODO: there must be a better way.
        fits = numpy.zeros(fitnesses.size(), dtype=numpy.double)
        cdef: 
            np.npy_double[:] np_fits = fits
            size_t i

        for i in range(fitnesses.size()):
            np_fits[i] = fitnesses[i]

        return fits

    property weighting:
        def __get__(self):
            return self._base.weighting

        def __set__(self, vector[double] wghts):
            self._base.set_weighting(wghts)

    property identifier:
        def __get__(self):
            return self._base.identifier

    property name:
        def __get__(self):
            return self._base.name

    def as_array(self):
        return numpy.array(self._base.optimal_rates)

    property scoring_method:
        def __get__(self):
            return self._base.scoring_method
        def __set__(self, ScoringMethod method):
            self._base.scoring_method = method

    property strength:
        def __get__(self):
            return self._base.strength
        def __set__(self, double s):
            self._base.strength = s

    def calc_categories(self):
        """Categorise the targets"""
        cat_dict = {}
        cats = []
        cat_n = 0
        for et in self.as_array():
            et = tuple(et)
            if et in cat_dict:
                cats.append(cat_dict[et])
            else:
                cat_dict[et] = cat_n
                cats.append(cat_n)
                cat_n += 1
        return cats
            
    def calc_distinct_outputs(self):
        out = set()
        for et in self.as_array():
            out.add(tuple(et))
        return out

    def _construct(self, init_func):
        a, b = self.world._this.cue_range

        # Slow and cumbersome, but it doesn't matter
        for i, e in enumerate(self.world.environments):
            # TODO: Clean up the refs here
            outputs = init_func(*e.as_array()[a:b])
            try:
                s = len(outputs)
            except TypeError:
                # Must be a single value...
                outputs = [outputs]
                s = 1

            if len(outputs) != self.world.out_channels:
                raise RuntimeError(
                    "return value of Target function must be length %s" \
                    % self.world.out_channels)

            for j, val in enumerate(outputs):
                self._base.optimal_rates[i][j] = float(val)


# NOTE: allow weighting to be None for backward compatibility
def _default_target(World w, name, ident, rates, weighting=None,
                      scoring_method=None, strength=None):
    t = DefaultTarget(w, None, name, ident=ident)
    # Manually construct these
    t._this.optimal_rates = rates
    if weighting is not None:
        t.weighting = weighting
    if scoring_method is not None:
        t.scoring_method = scoring_method
    if strength is not None:
        t.strength = strength
    return t

cdef class DefaultTarget(BaseTarget):
    def __cinit__(self, World w, init_func=None, name="", ident=-1,
                  scoring_method=0, strength=0.0):
        self.world = w
        self._this = new cDefaultTarget(w._shared, name, ident, 
                                        scoring_method, strength)
        self._base = self._this
        if init_func:
            self._construct(init_func)

    def __reduce__(self):
        return _default_target, (
            self.world, 
            self._this.name, 
            self._this.identifier,
            self._this.optimal_rates,
            self._this.weighting,
            self._this.scoring_method,
            self._this.strength)


def _noisy_target(World w, name, ident, rates, weighting,
                  scoring_method, strength, 
                  perturb_count, perturb_prop, env_only):
    t = NoisyTarget(w, None, name, ident, 
                    scoring_method, strength,
                    perturb_count, perturb_prop, env_only)
    # Manually construct these
    t._base.optimal_rates = rates
    t._base.weighting = weighting
    return t

            
cdef class NoisyTarget(BaseTarget):
    def __cinit__(self, World w, init_func=None, name="", ident=-1,
                  scoring_method=0, strength=0.0, 
                  perturb_count=1, perturb_prop=1.0, env_only=True):
        self.world = w
        self._this = new cNoisyTarget(
            w._shared, name, ident, 
            scoring_method, strength, perturb_count, perturb_prop, env_only)

        self._base = self._this
        if init_func:
            self._construct(init_func)

    def __reduce__(self):
        return _noisy_target, (
            self.world, 
            self._this.name, 
            self._this.identifier,
            self._this.optimal_rates,
            self._this.weighting,
            self._this.scoring_method,
            self._this.strength,
            self._this.perturb_count,
            self._this.perturb_prop,
            self._this.env_only,
        )

    property env_only:
        def __get__(self):
            return self._this.env_only

    property perturb_count:
        def __get__(self):
            return self._this.perturb_count

    property perturb_prop:
        def __get__(self):
            return self._this.perturb_prop


# NOTE: allow weighting to be None for backward compatibility
def _multi_target(World w, name, ident, rates, weighting=None,
                      scoring_method=None, strength=None):
    t = MultiTarget(w, None, name, ident=ident)
    # Manually construct these
    t._this.optimal_rates = rates
    if weighting is not None:
        t.weighting = weighting
    if scoring_method is not None:
        t.scoring_method = scoring_method
    if strength is not None:
        t.strength = strength
    return t


cdef class MultiTarget(BaseTarget):
    def __cinit__(self, World w, init_func=None, name="", ident=-1,
                  scoring_method=0, strength=0.0):
        self.world = w
        self._this = new cMultiTarget(w._shared, name, ident, 
                                        scoring_method, strength)
        self._base = self._this
        if init_func:
            self._construct(init_func)

        c = Channels(w)
        c.unchecked_set(2)
        self._this.pulses.push_back(c._this);

    def __reduce__(self):
        return _multi_target, (
            self.world, 
            self._this.name, 
            self._this.identifier,
            self._this.optimal_rates,
            self._this.weighting,
            self._this.scoring_method,
            self._this.strength)
