from utility cimport *
from grn cimport *

cdef class ChannelStateFrozen:
    cdef:
        cChannelState cchannel_state
        cFactory_ptr cfactory_ptr
    cdef init(self, cFactory_ptr &f, cChannelState &p)

cdef class ChannelState(ChannelStateFrozen):
    pass

cdef class Factory:
    cdef:
        cFactory_ptr cfactory_ptr
        cFactory *cfactory
        cGeneFactory *cgenefactory
        readonly:
            object params
            object cis_class
        object _environments

cdef class Target:
    cdef:
        cTarget *ctarget
        readonly:
            Factory factory

cdef class Network:
    cdef:
        readonly:
            Factory factory
            bint ready
            bint dirty

        # Because we hold a reference to the shared_ptr, we know we can always
        # safely access the ACTUAL pointer. We keep the pointer around too, as
        # it makes our life easier. The cost is a tiny bit of space.
        cNetwork_ptr ptr
        cNetwork *cnetwork
        object _genes, _attractors, _rates

    cdef bind_to(self, cNetwork_ptr &ptr)

cdef class Gene:
    cdef:
        # Assumption: Networks CANNOT mess with genes number once a network
        # has been established (You must copy and mutate a network).
        cGene *cgene
        object _modules

        readonly:
            # By holding this ref, we ensure the pointer is always valid (pace
            # what I said above. DON'T mess with the genes!)
            Network network
            size_t gene_number

cdef class CisModule:
    cdef:
        cCisModule *ccismodule
    
        readonly:
            Gene gene

    cdef reset_network(self)

cdef class NetworkAnalysis:
    cdef:
        cNetworkAnalysis *canalysis
        readonly:
            Network network

cdef class NetworkCollection:
    cdef:
        readonly: 
            Factory factory
        cNetworkVector cnetworks

    cdef object get_at(self, size_t i)

