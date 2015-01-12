from . import core
from .operand import Operand
from .logic2_ext import Constructor
from .core_ext import World, Target, Population, NetworkAnalysis
__all__ = ["World", "Constructor", "Parameters", "Operand",
           "Population", "Target", "NetworkAnalysis"]

class Parameters(core.Parameters):
    def __init__(self, **kwargs):
        # Add new parameters here-->
        self.operands = [Operand.AND, Operand.OR, Operand.NOT_A_AND_B,
                         Operand.A_AND_NOT_B, Operand.FALSE]
        self.factory_class = Constructor
        self.target_class = Target
        # <---
        core.Parameters.__init__(self, **kwargs)

