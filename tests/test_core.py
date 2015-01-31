import bricolage.core as T
import cPickle as pickle
import pathlib
import pytest

@pytest.fixture
def target_3x2():
    """Return a function for initialising a target that has 3 inputs and 2
    outputs"""
    def make_target(a, b, c):
        f1 = .5 if a and b or not c else 1.0
        f2 = 1 if ((a or c) and not (a and b)) and b else 0
        return f1, f2
    return make_target

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

def test_target(target_3x2):
    p = T.Parameters(
        cue_channels=3,
        reg_channels=3,
        out_channels=2,
    )
    w = T.World(p)
    t = T.Target(w, target_3x2)
    assert t.as_array().shape == (pow(2, 3), 2)


def test_pickling(tmpdir, target_3x2):
    tmpdir = pathlib.Path(str(tmpdir))
    p = T.Parameters(
        cue_channels=3,
        reg_channels=3,
        out_channels=2,
    )
    w = T.World(p)
    with open(str(tmpdir / 'world1.pickle'), 'wb') as f:
        pickle.dump(w, f, -1)

    with open(str(tmpdir / 'world1.pickle'), 'rb') as f:
        w2 = pickle.load(f)

    assert dir(w2.params) == dir(w.params)
    assert w.cue_channels == w2.cue_channels
    assert w.reg_channels == w2.reg_channels
    assert w.out_channels == w2.out_channels

    t = T.Target(w, target_3x2)
    with open(str(tmpdir / 'target1.pickle'), 'wb') as f:
        pickle.dump(t, f, -1)

    with open(str(tmpdir / 'target1.pickle'), 'rb') as f:
        t2 = pickle.load(f)

    assert (t.as_array() == t2.as_array()).all()
    assert t.name == t2.name
    assert t.identifier == t2.identifier


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

    # Now try with state setting
    ss = w.get_random_state()
    a = [w.get_random_double(0, 1) for _ in range(100)]

    w.set_random_state(ss)
    b = [w.get_random_double(0, 1) for _ in range(100)]

    assert a == b
        
    

