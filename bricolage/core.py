from core_ext import World

__all__ = ["World", "Parameters"]

class Parameters(object):
    def __init__(self, **kwargs):
        # Defaults are provided here
        self.seed = 1
        self.cis_count = 3
        self.gene_count = 3
        self.cue_channels = 2
        self.reg_channels = 1
        self.out_channels = 1

    def _override(self, kwargs):
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                raise RuntimeError("Invalid Setting: {}, in Parameters".format(k))
