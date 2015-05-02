import core

from .core_ext import World, Target, Population, SelectionModel
from .analysis_ext import NetworkAnalysis
from .threshold3_ext import Constructor

__all__ = ["World", "Constructor", "Parameters",
           "Population", "Target", "NetworkAnalysis", "SelectionModel"]

class Parameters(core.Parameters):
    def __init__(self, **kwargs):
        # Add new parameters here-->
        self.factory_class = Constructor
        self.target_class = Target
        # <---
        core.Parameters.__init__(self, **kwargs)

