import core
from .core_ext import (World, Network, Target, NetworkCollection,
                       NetworkAnalysis)
from .threshold3_ext import Constructor

__all__ = ["World", "Constructor", "Network", "Parameters",
           "NetworkCollection", "Target", "NetworkAnalysis"]

class Parameters(core.Parameters):
    def __init__(self, **kwargs):
        # Add new parameters here-->
        
        # <---
        core.Parameters.__init__(self, **kwargs)

