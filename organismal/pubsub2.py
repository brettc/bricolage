from pubsub2_ext import Factory, Products

# __all__ = 

class Parameters(object):
    def __init__(self, **kwargs):
        # Defaults are provided here
        self.gene_count = 10

        self.population_size = 20
        # self.cis_count = 3
        # self.cue_shapes = 2
        # self.reg_shapes = 3
        # self.out_shapes = 1
        # self.ops = 'and or not_x not_y'.split()
        # self.mutation_rate = .01

        self._override(kwargs)

        # Below here are internal calculated values
        self._init()

    def _override(self, kwargs):
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                pass
                # log.warning("'%s' is not a valid setting", k)

    def _init(self):
        self.blarg = None
        self.boodle = 10
        # self._calc_sizes()
        # self._init_ops()
        # self._init_envs()
        pass

