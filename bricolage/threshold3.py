import core
from core import MutateType
from .core_ext import World, Population
from .analysis_ext import NetworkAnalysis
from .targets_ext import DefaultTarget
from .core_ext import SelectionModel
from .threshold3_ext import Factory

__all__ = ["World", "Factory", "Parameters", "MutateType",
           "Population", "DefaultTarget", "NetworkAnalysis", "SelectionModel"]

class Parameters(core.Parameters):
    def __init__(self, **kwargs):
        # Add new parameters here-->
        self.factory_class = Factory
        self.target_class = DefaultTarget
        # <---
        core.Parameters.__init__(self, **kwargs)

