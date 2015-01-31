from .core_ext import SelectionModel, World, Target

__all__ = ["World", "SelectionModel", "Parameters", "Target"]

class Parameters(object):
    def __init__(self, **kwargs):
        # Defaults are provided here
        self.seed = 1
        self.cis_count = 3
        self.gene_count = 3
        self.cue_channels = 2
        self.reg_channels = 1
        self.out_channels = 1
        self.selection_class = SelectionModel
        self.population_size = 100

        self._override(kwargs)

    def _override(self, kwargs):
        # TODO: Something clever here
        for k, v in kwargs.items():
            setattr(self, k, v)
            # if hasattr(self, k):
    #         else:
    #             # TODO: Issue a warning
    #             pass
                # raise RuntimeError("Invalid Setting: {}, in Parameters".format(k))
