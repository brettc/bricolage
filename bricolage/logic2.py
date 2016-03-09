from . import core
from .operand import Operand
from .logic2_ext import Factory
from .core_ext import World, DefaultTarget, Population
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

