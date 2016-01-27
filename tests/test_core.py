import bricolage.core as T
import cPickle as pickle
import pathlib
import numpy


def make_target1(a, b, c):
    f1 = .5 if a and b or not c else 1.0
    f2 = 1 if ((a or c) and not (a and b)) and b else 0
    return f1, f2


def make_target2(a, b, c):
    f1 = .25 if (a or b) and (not a and not c) else 1.0
    f2 = 1 if ((a or b) and not (a and b)) and c else 0
    return f1, f2


def bowtie_target(a, b, c):
    if (a and not c) or (b and c):
        return [1, .5, .25]
    return [0, 0, 0]


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


def test_target():
    p = T.Parameters(
        cue_channels=3,
        reg_channels=3,
        out_channels=2,
    )
    w = T.World(p)
    t = T.Target(w, make_target1)
    assert t.as_array().shape == (pow(2, 3), 2)

    # Default
    assert t.weighting == [.5, .5]

    t.weighting = [1, 4]
    assert t.weighting == [.2, .8]


def test_pickling_world(tmpdir):
    tmpdir = pathlib.Path(str(tmpdir))
    p = T.Parameters(
        seed=99,
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
    assert w.get_random_state() == w2.get_random_state()
    assert w.next_network_id == w2.next_network_id
    assert w.next_target_id == w2.next_target_id


def test_pickling_target(tmpdir):
    tmpdir = pathlib.Path(str(tmpdir))
    p = T.Parameters(
        cue_channels=3,
        reg_channels=3,
        out_channels=2,
    )
    w = T.World(p)

    # Now ensure that pickling Targets works too
    t1 = T.Target(w, make_target1, name='a')
    assert t1.scoring_method == T.ScoringMethod.LINEAR
    assert t1.strength == 0.0

    t2 = T.Target(w, make_target2, name='b',
                  scoring_method=T.ScoringMethod.EXPONENTIAL, strength=4.0)
    t2.weighting = [1, 2]

    with open(str(tmpdir / 'target1.pickle'), 'wb') as f:
        pickle.dump((t1, t2), f, -1)

    with open(str(tmpdir / 'target1.pickle'), 'rb') as f:
        rt1, rt2 = pickle.load(f)

    assert (t1.as_array() == rt1.as_array()).all()
    assert (t2.as_array() == rt2.as_array()).all()

    assert t1.name == rt1.name
    assert t2.name == rt2.name

    assert t1.identifier == rt1.identifier
    assert t2.identifier == rt2.identifier

    assert t1.weighting == rt1.weighting
    assert t2.weighting == rt2.weighting

    assert t1.scoring_method == rt1.scoring_method
    assert t2.scoring_method == rt2.scoring_method

    assert t1.strength == rt1.strength
    assert t2.strength == rt2.strength


def test_scoring_methods(bowtie_database):
    pop = bowtie_database.population
    # Use different identifiers to force recalculation
    targ1 = T.Target(pop.factory.world, bowtie_target, ident=2)
    targ2 = T.Target(pop.factory.world, bowtie_target, ident=3,
                     scoring_method=T.ScoringMethod.EXPONENTIAL,
                     strength=1)
    targ3 = T.Target(pop.factory.world, bowtie_target, ident=4,
                     scoring_method=T.ScoringMethod.EXPONENTIAL_VEC,
                     strength=1)

    targ1.assess_collection(pop)
    f1 = pop.fitnesses
    targ2.assess_collection(pop)
    f2 = pop.fitnesses
    targ3.assess_collection(pop)
    f3 = pop.fitnesses
    ones1 = numpy.where(f1 == 1.0)[0]
    ones2 = numpy.where(f2 == 1.0)[0]
    ones3 = numpy.where(f3 == 1.0)[0]
    assert (ones1 == ones2).all()
    assert (ones1 == ones3).all()


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
