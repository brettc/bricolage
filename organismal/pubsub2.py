from operand import Operand
from pubsub2_ext import Factory, NetworkCollection

class Parameters(object):
    def __init__(self, **kwargs):
        # Defaults are provided here
        self.seed = 1
        self.population_size = 20
        self.cis_count = 3
        self.cue_channels = 2
        self.reg_channels = 3
        self.out_channels = 1
        self.mutation_rate = .01
        self.operands = [Operand.AND, Operand.OR, Operand.A_AND_NOT_B]

        self._override(kwargs)

        # Below here are internal calculated values
        self._init()
        self._validate()

    def _override(self, kwargs):
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                raise RuntimeError("Invalid Operations")

    def _init(self):
        # Calculate the total number of elements given the overlap
        # An example, with an overlap of 2, looks like this:
        # cue_range [ 0 1 2 3 ]
        # sub_range [ 0 1 2 3 4 5 ]
        # pub_range         [ 4 5 6 7 8 9 10 ]
        # out_range             [ 6 7 8 9 10 ]
        # reg_range         [ 4 5 ] 
        
        self.reg_gene_count = self.reg_channels
        self.sub_channels = self.cue_channels + self.reg_channels

        self.channel_count = self.cue_channels + self.reg_channels + \
            self.out_channels

        self.sub_range = 0, self.sub_channels
        self.reg_range = self.cue_channels, self.sub_channels
        self.pub_range = self.cue_channels, self.sub_channels + self.out_channels
        self.cue_range = 0, self.cue_channels
        self.out_range = self.sub_channels, self.channel_count

        # Let's make these, they'll be useful
        self.cue_signals = range(*self.cue_range)
        self.reg_signals = range(*self.reg_range)
        self.out_signals = range(*self.out_range)

        self.sub_signals = range(*self.sub_range)
        self.pub_signals = range(*self.pub_range)

        # Total gene count requires 
        self.gene_count = self.reg_gene_count + self.out_channels
        self.env_count = pow(2, self.cue_channels)

        # This is per CIS module, rather than per site
        self.individual_mutation_rate = self.gene_count * self.cis_count * self.mutation_rate

        # Lastly, work out our char_size
        # Should never need more than 2, as 99 is really big
        if len(self.cue_signals) <= 10 and \
                len(self.reg_signals) <= 10 and \
                len(self.out_signals) <= 10:
            self.char_size = 2
        else:
            self.char_size = 3

    def name_for_channel(self, s):
        sz = self.char_size - 1
        # We use "mathematical" number, starting at 1.0
        if s >= self.out_range[0]:
            return "P{0:0{1:}d}".format(s + 1 - self.out_range[0], sz)

        if s >= self.reg_range[0]:
            return "T{0:0{1:}d}".format(s + 1 - self.reg_range[0], sz)

        return "E{0:0{1:}d}".format(s + 1 - self.cue_range[0], sz)

    def _validate(self):
        pass
