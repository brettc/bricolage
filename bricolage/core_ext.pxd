from utility cimport *
from core cimport *

cdef class ChannelStateFrozen:
    cdef:
        cChannelState _this
        cWorld_ptr world
    cdef init(self, cWorld_ptr &f, cChannelState &p)

cdef class ChannelState(ChannelStateFrozen):
    pass

cdef class World:
    cdef:
        cWorld_ptr _shared
        cWorld *_this
        object _params
        readonly:
            object sub_signals, pub_signals
            object cue_signals, reg_signals, out_signals
            object reserved_signals
        object _environments

cdef class Factory:
    cdef:
        cFactory_ptr _shared
        cFactory *_this
        int _secret_key
        readonly:
            World world
            object gene_class, module_class, network_class

cdef class BaseTarget:
    cdef:
        cBaseTarget *_base
        readonly:
            World world

cdef class DefaultTarget(BaseTarget):
    cdef:
        cDefaultTarget *_this

cdef class NoisyTarget(BaseTarget):
    cdef:
        cNoisyTarget *_this

cdef class Network:
    cdef:
        readonly:
            Factory factory

        # Because we hold a reference to the shared_ptr, we know we can always
        # safely access the ACTUAL pointer. We keep the pointer around too, as
        # it makes our life easier. The cost is a tiny bit of space.
        cNetwork_ptr _shared
        cNetwork *_this
        object _genes, _attractors, _rates

    cdef bind_to(self, cNetwork_ptr ptr)

    cdef _make_python_attractor(self, cChannelStateVector &attr)
    cdef _make_python_attractors(self, cAttractors &attrs)
    cdef _make_python_rates(self, cRatesVector &)
    cdef _make_python_rate(self, cRates &rates)

cdef class Gene:
    cdef:
        # Assumption: Networks CANNOT mess with genes number once a network
        # has been established (You must copy and mutate a network).
        # TODO: ensure this using point to const!
        cGene *_this
        object _modules

        readonly:
            # By holding this ref, we ensure the pointer is always valid (pace
            # what I said above. DON'T mess with the genes!)
            Network network
            size_t gene_number

cdef class CisModule:
    cdef:
        cCisModule *_this
        readonly:
            Gene gene

cdef class SelectionModel:
    cdef:
        cSelectionModel *_this
        readonly:
            World world

cdef class CollectionBase:
    cdef:
        cNetworkVector *_collection
        readonly:
            Factory factory

    cdef object get_at(self, size_t i)

cdef class Collection(CollectionBase):
    cdef:
        cNetworkVector *_this

cdef class Ancestry(CollectionBase):
    cdef:
        cNetworkVector *_this

cdef class Population(CollectionBase):
    cdef:
        cPopulation *_this

