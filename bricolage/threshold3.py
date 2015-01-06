import core
from .core_ext import (World, Target, Population,
                       NetworkAnalysis, SelectionModel)
from .threshold3_ext import Constructor

__all__ = ["World", "Constructor", "Parameters",
           "Population", "Target", "NetworkAnalysis", "SelectionModel"]

class Parameters(core.Parameters):
    def __init__(self, **kwargs):
        # Add new parameters here-->
        
        # <---
        core.Parameters.__init__(self, **kwargs)

