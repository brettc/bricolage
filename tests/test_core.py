# import pytest
from organismal import core_ext

def test_channels():
    c = core_ext.ChannelDef(2, 2, 2)
    print c.cue_channels
    print c.cue_range
