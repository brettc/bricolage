from .core_ext import NetworkCollection, NetworkAnalysis, Target, World, Network
from .threshold3_ext import Constructor

__all__ = [
    "Target", "NetworkCollection", "NetworkAnalysis", "World", "Network",
    "Parameters", "Constructor"
]

class Parameters(object):
    def __init__(self, **kwargs):
        # Defaults are provided here
        self.seed = 1
        self.cis_count = 3
        self.cue_channels = 2
        self.reg_channels = 3
        self.out_channels = 1

        self._override(kwargs)

    def _override(self, kwargs):
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                raise RuntimeError("Invalid Operations")

