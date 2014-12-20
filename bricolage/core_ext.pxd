from utility cimport *
from grn cimport *

cdef class ChannelStateFrozen:
    cdef:
        cChannelState cchannel_state
        cWorld_ptr cworld_ptr
    cdef init(self, cWorld_ptr &f, cChannelState &p)

cdef class ChannelState(ChannelStateFrozen):
    pass

cdef class World:
    cdef:
        cWorld_ptr cworld_ptr
        cWorld *cworld
        readonly:
            object sub_signals, pub_signals
            object cue_signals, reg_signals, out_signals
            object reserved_signals
        object gene_class, module_class
        object _environments

# cdef class Target:
#     cdef:
#         cTarget *ctarget
#         readonly:
#             World world

cdef class Network:
    cdef:
        readonly:
            World world

        # Because we hold a reference to the shared_ptr, we know we can always
        # safely access the ACTUAL pointer. We keep the pointer around too, as
        # it makes our life easier. The cost is a tiny bit of space.
        cNetwork_ptr ptr
        cNetwork *cnetwork
        object _genes, _attractors, _rates
        object gene_class

    cdef bind_to(self, cNetwork_ptr &ptr)

cdef class Gene:
    cdef:
        # Assumption: Networks CANNOT mess with genes number once a network
        # has been established (You must copy and mutate a network).
        # TODO: ensure this using point to const!
        cGene *cgene
        object _modules

        readonly:
            # By holding this ref, we ensure the pointer is always valid (pace
            # what I said above. DON'T mess with the genes!)
            Network network
            size_t gene_number
#
cdef class CisModule:
    cdef:
        cCisModule *ccismodule
        readonly:
            Gene gene

cdef class NetworkCollection:
    cdef:
        readonly: 
            World world
        cNetworkVector cnetworks

    cdef object get_at(self, size_t i)

#
# cdef class NetworkAnalysis:
#     cdef:
#         cNetworkAnalysis *canalysis
#         readonly:
#             Network network
#
