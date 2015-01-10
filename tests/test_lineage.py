import pytest
import bricolage.logic2 as T
import bricolage.lineage as L
import pathlib

@pytest.fixture
def p_3x2():
    o = T.Operand
    ops = o.NOT_A_AND_B, o.A_AND_NOT_B, o.NOR, o.AND
    return T.Parameters(seed=1, operands=ops, cis_count=2, reg_channels=5,
                        out_channels=2, cue_channels=3, population_size=1000)
@pytest.fixture
def c_3x2(p_3x2):
    world = T.World(p_3x2)
    return T.Constructor(world, p_3x2)

@pytest.fixture
def target_3x2():
    """Return a function for initialising a target that has 3 inputs and 2
    outputs"""
    def make_target(a, b, c):
        f1 = .5 if a and b or not c else 1.0
        f2 = 1 if ((a or c) and not (a and b)) and b else 0
        return f1, f2
    return make_target

def test_numpy_export(c_3x2):
    p1 = T.Population(c_3x2, 1000)
    as_array = c_3x2.to_numpy(p1)

    assert as_array.dtype == c_3x2.dtype()

    p2 = T.Population(c_3x2, 0)
    c_3x2.from_numpy(as_array, p2)

    # Did we get exactly the same thing back?
    assert p1.size == p2.size
    for n1, n2 in zip(p1, p2):
        assert n1.attractors == n2.attractors
        assert (n1.rates == n2.rates).all()
        for g1, g2 in zip(n1.genes, n2.genes):
            assert g1.pub == g2.pub
            for m1, m2 in zip(g1.modules, g2.modules):
                assert m1.op == m2.op
                assert m1.channels == m2.channels


def test_creation(tmpdir):
    base = pathlib.Path(str(tmpdir))
    path = base / 'creation.db'
    params = T.Parameters(population_size=10)
    a = L.Lineage(path, params, T.Constructor)
    del a
    b = L.Lineage(path)
    del b


def assert_pops_equal(p1, p2):
    assert p1.size == p2.size
    for n1, n2 in zip(p1, p2):
        assert n1.identifier == n2.identifier
        assert n1.parent_identifier == n2.parent_identifier
        assert n1.generation == n2.generation
        assert n1.attractors == n2.attractors
        assert (n1.rates == n2.rates).all()
        for g1, g2 in zip(n1.genes, n2.genes):
            assert g1.pub == g2.pub
            for m1, m2 in zip(g1.modules, g2.modules):
                assert m1.op == m2.op
                assert m1.channels == m2.channels

def test_repeatability(p_3x2, target_3x2):
    """Make sure that the same seed produces the same outcome"""
    a = L.SnapshotLineage(params=p_3x2, factory_class=T.Constructor)
    target = T.Target(a.world, target_3x2)
    sel = T.SelectionModel(a.world)
    for i in range(100):
        a.assess(target)
        a.next_generation(.001, sel)

    b = L.SnapshotLineage(params=p_3x2, factory_class=T.Constructor)
    target = T.Target(b.world, target_3x2)
    sel = T.SelectionModel(b.world)
    for i in range(100):
        b.assess(target)
        b.next_generation(.001, sel)

    assert_pops_equal(a.population, b.population)

def test_snapshot_reload(p_3x2, target_3x2, tmpdir):
    base = pathlib.Path(str(tmpdir))
    path = base / 'reload.db'
    a = L.SnapshotLineage(params=p_3x2, factory_class=T.Constructor)
    a.save_snapshot(path)
    p1 = a.population
    del a
    # Now reload and compare the population
    b = L.SnapshotLineage(path=path)
    assert_pops_equal(p1, b.population)
    del b

def test_snapshot_lineage(p_3x2, target_3x2, tmpdir):
    base = pathlib.Path(str(tmpdir))
    times = 5
    generations = 20
    mrate = .001

    # Generate a bunch of snapshots
    a = L.SnapshotLineage(params=p_3x2, factory_class=T.Constructor)
    target = T.Target(a.world, target_3x2)
    sel = T.SelectionModel(a.world)
    for i in range(times):
        for j in range(generations):
            a.assess(target)
            a.next_generation(mrate, sel)
        path = base / 'snapshot-{}.db'.format(i)
        a.save_snapshot(path)
    del a

    # Same seed will generate the same lineage
    b = L.SnapshotLineage(params=p_3x2, factory_class=T.Constructor)
    target = T.Target(b.world, target_3x2)
    sel = T.SelectionModel(b.world)
    for i in range(times):
        for j in range(generations):
            b.assess(target)
            b.next_generation(mrate, sel)
        path = base / 'snapshot-{}.db'.format(i)
        c = L.SnapshotLineage(path=path)

        # Reloads should be the same
        assert_pops_equal(b.population, c.population)

def test_restarting(tmpdir, p_3x2, target_3x2):
    fname = pathlib.Path(str(tmpdir)) / 'selection.db'

    # Set it all up
    a = L.Lineage(fname, p_3x2, T.Constructor)
    sel = T.SelectionModel(a.world)
    target = T.Target(a.world, target_3x2)

    # Run selection for 100 generations
    for i in range(100):
        a.assess(target)
        a.next_generation(.01, sel)
    x = [net.identifier for net in a.population]

    # Kill off the reference (automatically closing the file)
    del a

    # Reload it and delete, but save the idents
    b = L.Lineage(fname)
    y = [net.identifier for net in b.population]
    del b

    # It should be the same
    assert x == y

    # Reload again, this time running selection further.
    # Save the 100th generation too.
    c = L.Lineage(fname)
    p1 = c.get_generation(100)
    for i in range(100):
        c.assess(target)
        c.next_generation(.01, sel)
    del c

    # Relead yet again. Reload the 100th generation again. Check it is exactly
    # the same
    d = L.Lineage(fname)
    p2 = d.get_generation(100)
    z = [net.identifier for net in p2]
    assert y == z

    # This means everything is the same!
    assert_pops_equal(p1, p2)

    # Pull out an ancestry
    anc = d.get_ancestry(z[0])
    
    # NOTE: should pull this out separately
    prev_i = -1
    for n in anc:
        i = n.identifier
        p = n.parent_identifier
        assert p == prev_i
        prev_i = i

# def test_ancestry():
#     assert 1 == 2
    

