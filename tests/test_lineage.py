import pytest
import bricolage.logic2 as T
import bricolage.lineage as L
import pathlib

@pytest.fixture
def p_3x2():
    o = T.Operand
    ops = o.NOT_A_AND_B, o.A_AND_NOT_B, o.NOR, o.AND
    return T.Parameters(seed=1, operands=ops, cis_count=2, reg_channels=5,
                        out_channels=2, cue_channels=3, population_size=100)
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
    path = base / 'create.db'
    params = T.Parameters(population_size=10)
    a = L.FullLineage(path, params, T.Constructor)
    del a
    # b = L.Lineage(path)
    # del b


def assert_pops_equal(p1, p2):
    assert p1.size == p2.size
    for n1, n2 in zip(p1, p2):
        assert n1.identifier == n2.identifier
        assert n1.parent_identifier == n2.parent_identifier
        assert n1.generation == n2.generation
        for g1, g2 in zip(n1.genes, n2.genes):
            assert g1.pub == g2.pub
            for m1, m2 in zip(g1.modules, g2.modules):
                assert m1.op == m2.op
                assert m1.channels == m2.channels
        # No blank attractors allowed
        assert n1.attractors
        assert n1.attractors == n2.attractors
        assert (n1.rates == n2.rates).all()

def test_repeatability(tmpdir, p_3x2, target_3x2):
    """Make sure that the same seed produces the same outcome"""
    base = pathlib.Path(str(tmpdir))
    path = base / 'r1.db'
    a = L.SnapshotLineage(path, params=p_3x2, factory_class=T.Constructor)
    target = T.Target(a.world, target_3x2)
    sel = T.SelectionModel(a.world)
    for i in range(100):
        a.assess(target)
        a.next_generation(.001, sel)

    path = base / 'r2.db'
    b = L.SnapshotLineage(path, params=p_3x2, factory_class=T.Constructor)
    target = T.Target(b.world, target_3x2)
    sel = T.SelectionModel(b.world)
    for i in range(100):
        b.assess(target)
        b.next_generation(.001, sel)

    assert_pops_equal(a.population, b.population)
    del a
    del b

def test_snapshot_reload(p_3x2, target_3x2, tmpdir):
    base = pathlib.Path(str(tmpdir))
    path = base / 'reload.db'
    a = L.SnapshotLineage(path, params=p_3x2, factory_class=T.Constructor)
    a.save_snapshot()
    p1 = a.population
    r1 = a.world.get_random_state()
    del a
    # Now reload and compare the population
    b = L.SnapshotLineage(path=path)
    assert_pops_equal(p1, b.population)
    assert b.world.get_random_state() == r1
    del b

def test_snapshot_lineage(p_3x2, target_3x2, tmpdir):
    base = pathlib.Path(str(tmpdir))
    path_1 = base / 'snap1.db'
    times = 5
    generations = 20
    mrate = .001

    # Generate a bunch of snapshots
    a = L.SnapshotLineage(path_1, params=p_3x2, factory_class=T.Constructor)
    a_target = T.Target(a.world, target_3x2)
    a_sel = T.SelectionModel(a.world)
    for i in range(times):
        for j in range(generations):
            a.assess(a_target)
            a.next_generation(mrate, a_sel)
        a.save_snapshot()

    # Delete and reload a
    del a
    a = L.SnapshotLineage(path_1) 
    a_target = T.Target(a.world, target_3x2)
    a_sel = T.SelectionModel(a.world)

    # Same seed will generate the same lineage
    path_2 = base / 'snap2.db'
    b = L.SnapshotLineage(path_2, params=p_3x2, factory_class=T.Constructor)
    b_target = T.Target(b.world, target_3x2)
    b_sel = T.SelectionModel(b.world)
    for i in range(times):
        for j in range(generations):
            b.assess(b_target)
            b.next_generation(mrate, b_sel)
        p1 = a.get_generation(b.generation)

        # Reloads should be the same
        assert_pops_equal(b.population, p1)

    # Now run them in parallel -- restarting should be the same as if we just
    # ran it from the beginning
    assert_pops_equal(a.population, b.population)
    assert a.world.get_random_state() == b.world.get_random_state()

    for i in range(20):
        b.assess(b_target)
        b.next_generation(mrate, b_sel)
        a.assess(a_target)
        a.next_generation(mrate, a_sel)

    assert_pops_equal(a.population, b.population)
    assert a.world.get_random_state() == b.world.get_random_state()

def test_full_lineage(tmpdir, p_3x2, target_3x2):
    path = pathlib.Path(str(tmpdir)) / 'selection.db'

    # Set it all up
    a = L.FullLineage(path, p_3x2, T.Constructor)
    sel = T.SelectionModel(a.world)
    target = T.Target(a.world, target_3x2)

    # Run selection for 100 generations
    for i in range(100):
        a.assess(target)
        a.next_generation(.01, sel)
    pa = a.population

    # Kill off the reference (automatically closing the file)
    del a

    # Reload it and delete, but save the idents
    b = L.FullLineage(path)
    pb = b.population
    assert_pops_equal(pa, pb)
    # print b

    # Reload again, this time running selection further.
    # Save the 100th generation too.
    c = L.FullLineage(path)
    sel = T.SelectionModel(c.world)
    target = T.Target(c.world, target_3x2)
    p1 = c.get_generation(100)
    for i in range(100):
        c.assess(target)
        c.next_generation(.01, sel)
    del c

    # Relead yet again. Reload the 100th generation again. Check it is exactly
    # the same
    d = L.FullLineage(path)
    p2 = d.get_generation(100)
    assert_pops_equal(p1, p2)

    # This means everything is the same!
    # assert_pops_equal(p1, p2)
    # anc = b.get_ancestry(pb[0].identifier)
    # print anc

# def test_ancestry():
    # Pull out an ancestry
    # anc = d.get_ancestry(z[0])
    
    # NOTE: should pull this out separately
    # prev_i = -1
    # for n in anc:
    #     i = n.identifier
    #     p = n.parent_identifier
    #     assert p == prev_i
    #     prev_i = i
    

