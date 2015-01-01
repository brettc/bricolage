# import pytest
import bricolage.core as T

def test_world():
    cue = 3
    reg = 4
    out = 3
    p = T.Parameters(
        cue_channels=cue,
        reg_channels=reg,
        out_channels=out,
    )
    w = T.World(p)
    assert w.cue_channels == cue
    assert w.reg_channels == reg
    assert w.out_channels == out
    assert w.channel_count == 2 + cue + reg + out

def test_channelstate():
    p = T.Parameters(
        cue_channels=3,
        reg_channels=4,
        out_channels=3,
    )
    w = T.World(p)

    e2 = w.environments[-1]
    e2_again = w.environments[-1]

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

def test_random_engine():
    p = T.Parameters(
        cue_channels=3,
        reg_channels=4,
        out_channels=3,
    )
    w = T.World(p)
    w.seed_random_engine(1)
    first_time = [w.get_random_double(0, 1) for _ in range(20)]
    first_time += [w.get_random_int(0, 100) for _ in range(20)]
    
    w.seed_random_engine(1)
    second_time = [w.get_random_double(0, 1) for _ in range(20)]
    second_time += [w.get_random_int(0, 100) for _ in range(20)]

    assert first_time == second_time
