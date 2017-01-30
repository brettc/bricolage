from . import core
import numpy as np
from .operand import Operand
from .logic2_ext import Factory
from .core_ext import World, Population
from .targets_ext import DefaultTarget
from bricolage.core_ext import Collection
__all__ = ["World", "Factory", "Parameters", "Operand",
           "Population", "DefaultTarget"] 

class Parameters(core.Parameters):
    def __init__(self, **kwargs):
        # Add new parameters here-->
        self.operands = [Operand.AND, Operand.OR, Operand.NOT_A_AND_B,
                         Operand.A_AND_NOT_B, Operand.FALSE]
        self.factory_class = Factory
        self.target_class = DefaultTarget
        # <---
        core.Parameters.__init__(self, **kwargs)


def count_diff_networks(net1, net2):
    f1 = net1.factory
    assert net2.factory is f1
    coll = Collection(f1)
    coll.add(net1)
    coll.add(net2)

    arr = f1.to_numpy(coll)
    subs = arr['sub']
    sub_diffs = subs[0] != subs[1]

    # Check to see if there ops changes where the subs did not change
    binding_diffs = sub_diffs.sum(axis=2) != 0
    mods_diff = np.where(binding_diffs != 0)
    print mods_diff

    ops = arr['op']
    op_diffs = ops[0] != ops[1]
    #
    # count of CRMs where binding sites are the same but the logic has
    # changed.
    # ops_only_count = (op_diffs & binding_diffs).sum()

    # Number of individual changes to the binding sites (we assume the logic
    # has changed)
    # op_count = op_diffs.sum()
    # sub_count = sub_diffs.sum()
    #
    # return sub_count, op_count


def modules_changed(net1, net2):
    f1 = net1.factory
    assert net2.factory is f1
    coll = Collection(f1)
    coll.add(net1)
    coll.add(net2)

    arr = f1.to_numpy(coll)
    subs = arr['sub']
    sub_diffs = subs[0] != subs[1]

    # Check to see if there ops changes where the subs did not change
    binding_diffs = sub_diffs.sum(axis=2) != 0

    ops = arr['op']
    ops_diffs = ops[0] != ops[1]

    any_mod_diffs = binding_diffs | ops_diffs

    return zip(*np.where(any_mod_diffs))

