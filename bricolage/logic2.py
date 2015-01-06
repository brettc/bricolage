from .operand import Operand
import core
from .core_ext import (World, Target, Population,
                       NetworkAnalysis, SelectionModel)
from .logic2_ext import Constructor

__all__ = ["World", "Constructor", "Parameters", "Operand",
           "Population", "Target", "NetworkAnalysis"]

class Parameters(core.Parameters):
    def __init__(self, **kwargs):
        # Add new parameters here-->
        self.operands = [Operand.AND, Operand.OR]
        # <---
        core.Parameters.__init__(self, **kwargs)

