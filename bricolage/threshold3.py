import core
from core import MutateType

from .core_ext import World, Target, Population, SelectionModel
from .analysis_ext import NetworkAnalysis
from .threshold3_ext import Factory

__all__ = ["World", "Factory", "Parameters", "MutateType", 
           "Population", "Target", "NetworkAnalysis", "SelectionModel"]

class Parameters(core.Parameters):
    def __init__(self, **kwargs):
        # Add new parameters here-->
        self.factory_class = Factory
        self.target_class = Target
        # <---
        core.Parameters.__init__(self, **kwargs)

