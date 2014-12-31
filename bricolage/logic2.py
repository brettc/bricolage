from .operand import Operand
import core
from core_ext import World, Network, NetworkCollection, Target, NetworkAnalysis
from logic2_ext import Constructor

__all__ = ["World", "Constructor", "Network", "Parameters", "Operand",
           "NetworkCollection", "Target", "NetworkAnalysis"]

class Parameters(core.Parameters):
    def __init__(self, **kwargs):
        core.Parameters.__init__(self, **kwargs)
        self.operands = [Operand.AND, Operand.OR]
        self._override(kwargs)

