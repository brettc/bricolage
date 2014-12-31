# import pytest
from bricolage.core import *

def test_world():
    w = World(Parameters())
    print w

def test_environments():
    # There should be X
    pass

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
        self.out_channels + 1

    self.cue_range = 1, self.cue_channels + 1
    self.sub_range = 1, self.sub_channels + 1
    self.reg_range = self.cue_channels + 1, self.sub_channels + 1
    self.pub_range = self.cue_channels + 1, self.sub_channels + \
        self.out_channels + 1
    self.out_range = self.sub_channels + 1, self.channel_count + 1

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
    # self.individual_mutation_rate = self.gene_count * self.cis_count * self.mutation_rate

    # Lastly, work out our char_size
    # Should never need more than 2, as 99 is really big
    if len(self.cue_signals) <= 10 and \
            len(self.reg_signals) <= 10 and \
            len(self.out_signals) <= 10:
        self.char_size = 2
    else:
        self.char_size = 3

def test_channelstate(p_3x2):
    f = T.Factory(p_3x2)
    e2 = f.environments[-1]
    e2_again = f.environments[-1]

    # We should get the same channels states out.
    assert e2 == e2_again
    assert e2 is e2_again

    # When we copy, they should be the same, but not identical.
    copy_e2 = e2.copy()
    assert e2 == copy_e2
    assert e2 is not copy_e2

    # Modify the state -- testing still work
    copy_e2.flip(0)
    assert e2 != copy_e2
    copy_e2.flip(0)
    assert e2 == copy_e2
