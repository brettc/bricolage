import pytest
from bricolage.operand import Operand
from bricolage import logic2
from bricolage import threshold3
import bricolage.lineage as L
import pathlib

@pytest.fixture(scope="module", params=[logic2, threshold3])
def module(request):
    return request.param

# Used for simple testing of one module
# @pytest.fixture
# def module():
#     return threshold3

@pytest.fixture
def p_3x2(module):
    o = Operand
    ops = o.NOT_A_AND_B, o.A_AND_NOT_B, o.NOR, o.AND
    return module.Parameters(
        seed=1, 
        operands=ops, 
        cis_count=2, 
        reg_channels=10,
        out_channels=2, 
        cue_channels=3, 
        population_size=100,
        mutation_rate=.001,
        replicates=10,
    )

def target1(a, b, c):
    f1 = .5 if a and b or not c else 1.0
    f2 = 1 if ((a or c) and not (a and b)) and b else 0
    return f1, f2

def target2(a, b, c):
    f1 = .5 if (a or b) and not c else 0
    f2 = 1 if ((a or c) and not (a and c)) or b else 0.5
    return f1, f2

def test_numpy_export(p_3x2, module):
    world = module.World(p_3x2)
    factory = module.Constructor(world)
    
    p1 = module.Population(factory, 1000)
    as_array = factory.to_numpy(p1)

    assert as_array.dtype == factory.dtype()

    p2 = module.Population(factory, 0)
    factory.from_numpy(as_array, p2)

    # Did we get exactly the same thing back?
    assert p1.size == p2.size
    for n1, n2 in zip(p1, p2):
        assert n1.attractors == n2.attractors
        assert (n1.rates == n2.rates).all()
        for g1, g2 in zip(n1.genes, n2.genes):
            assert g1.pub == g2.pub
            for m1, m2 in zip(g1.modules, g2.modules):
                assert m1.same_as(m2)
                assert m1.channels == m2.channels

def test_creation(tmpdir, module):
    base = pathlib.Path(str(tmpdir))
    path = base / 'create.db'
    params = module.Parameters(population_size=10)
    a = L.FullLineage(path, params)
    del a

def assert_pops_equal(p1, p2):
    assert p1.size == p2.size
    for n1, n2 in zip(p1, p2):
        assert n1.identifier == n2.identifier
        assert n1.parent_identifier == n2.parent_identifier
        assert n1.generation == n2.generation
        for g1, g2 in zip(n1.genes, n2.genes):
            assert g1.pub == g2.pub
            for m1, m2 in zip(g1.modules, g2.modules):
                assert m1.same_as(m2)

        # No blank attractors allowed
        assert n1.attractors
        assert n1.attractors == n2.attractors
        assert (n1.rates == n2.rates).all()

def test_snapshot_reload(p_3x2, tmpdir):
    base = pathlib.Path(str(tmpdir))
    path = base / 'reload.db'
    a = L.SnapshotLineage(path, params=p_3x2)
    a.add_target(target1)

    r1 = a.world.get_random_state()
    p1 = a.population

    a.save_snapshot()
    del a

    # Now reload and compare the population
    b = L.SnapshotLineage(path=path)
    assert_pops_equal(p1, b.population)
    r2 = b.world.get_random_state()
    assert r1 == r2
    assert b.world == b.targets[0].world
    assert b.world == b.factory.world
    del b

def test_snapshot_autosave(tmpdir, p_3x2):
    path = pathlib.Path(str(tmpdir)) / 'autosave.db'

    # Set it all up
    a = L.SnapshotLineage(path, p_3x2)
    a.add_target(target1)

    # Run selection for 100 generations
    for i in range(100):
        a.next_generation()
    pa = a.population
    # No explicit save!
    del a

    b = L.SnapshotLineage(path)
    assert_pops_equal(pa, b.population)
    b.close()

def test_snapshot_lineage(p_3x2, tmpdir):
    base = pathlib.Path(str(tmpdir))
    path_1 = base / 'snap1.db'
    times = 5
    generations = 20

    # Generate a bunch of snapshots
    a = L.SnapshotLineage(path_1, params=p_3x2)
    a.add_target(target1)
    for i in range(times):
        for j in range(generations):
            a.next_generation()
        a.save_snapshot()

    # Delete and reload a
    del a
    a = L.SnapshotLineage(path_1) 

    # Same seed will generate the same lineage
    path_2 = base / 'snap2.db'
    b = L.SnapshotLineage(path_2, params=p_3x2)
    b.add_target(target1)
    for i in range(times):
        for j in range(generations):
            b.next_generation()
        p1 = a.get_generation(b.generation)

        # Reloads should be the same
        assert_pops_equal(b.population, p1)

    # Now run them in parallel -- restarting should be the same as if we just
    # ran it from the beginning
    assert_pops_equal(a.population, b.population)
    assert a.world.get_random_state() == b.world.get_random_state()

    for i in range(20):
        b.next_generation()
        a.next_generation()

    assert_pops_equal(a.population, b.population)
    assert a.world.get_random_state() == b.world.get_random_state()

def test_repeatability(tmpdir, p_3x2):
    """Make sure that the same seed produces the same outcome"""
    base = pathlib.Path(str(tmpdir))
    path = base / 'r1.db'
    a = L.SnapshotLineage(path, params=p_3x2)
    a.add_target(target1)
    for i in range(100):
        a.next_generation()

    path = base / 'r2.db'
    b = L.SnapshotLineage(path, params=p_3x2)
    b.add_target(target1)
    for i in range(100):
        b.next_generation()

    assert_pops_equal(a.population, b.population)
    del a
    del b

def test_full_lineage(tmpdir, p_3x2):
    path = pathlib.Path(str(tmpdir)) / 'selection.db'

    # Set it all up
    a = L.FullLineage(path, p_3x2)
    a.add_target(target1)

    # Run selection for 100 generations
    for i in range(100):
        a.next_generation()
    pa = a.population

    # Kill off the reference (automatically closing the file)
    del a

    # Reload it and delete, but save the idents
    b = L.FullLineage(path)
    pb = b.population
    assert_pops_equal(pa, pb)

    # Reload again, this time running selection further.
    # Save the 100th generation too.
    c = L.FullLineage(path)
    c.add_target(target1)
    p1 = c.get_generation(100)
    for i in range(100):
        c.next_generation()
    del c

    # Relead yet again. Reload the 100th generation again. Check it is exactly
    # the same
    d = L.FullLineage(path)
    p2 = d.get_generation(100)
    assert_pops_equal(p1, p2)

    print {n.identifier: n for n in p2.get_best()}.items()

def test_ancestry(tmpdir, p_3x2):
    path = pathlib.Path(str(tmpdir)) / 'ancestry.db'

    # Set it all up
    a = L.FullLineage(path, p_3x2)
    a.add_target(target1)
    N = 1000
    for i in range(N):
        a.next_generation()

    del a
    b = L.FullLineage(path)

    # Pull out an ancestry
    anc = b.get_ancestry(b.population[0].identifier)

    # NOTE: should pull this out separately
    prev_i = -1
    for n in anc:
        i = n.identifier
        p = n.parent_identifier
        assert p == prev_i
        prev_i = i

def test_targets(p_3x2, tmpdir):
    base = pathlib.Path(str(tmpdir))
    path = base / 'targets.db'
    a = L.FullLineage(path, params=p_3x2)
    a.add_target(target1)
    for i in range(100):
        a.next_generation()
    
    pa1_gen = a.generation
    fa1 = [n.fitness for n in a.population]

    for i in range(100):
        a.next_generation()

    # Try something new
    a.add_target(target2)
    for i in range(100):
        a.next_generation()

    fa2 = [n.fitness for n in a.population]
    del a

    b = L.FullLineage(path=path)
    fb2 = [n.fitness for n in b.population]
    assert fa2 == fb2
    
    # This should automatically reload and recalculate the fitnesses using the
    # appropriate target
    fb1 = [n.fitness for n in b.get_generation(pa1_gen)]
    assert fa1 == fb1

def test_treatment(tmpdir, p_3x2):
    tmpdir = pathlib.Path(str(tmpdir))
    name = 'treat'
    path = tmpdir / name
    apath = tmpdir / (name + '_analysis')
    treat = L.Treatment(path, p_3x2, analysis_path=apath)

    def callback(rep, max_gen):
        lineage = rep.get_lineage()
        lineage.add_target(target1)
        while lineage.generation < max_gen:
            lineage.next_generation()
        assert rep.analysis_path is not rep.path
        assert rep.analysis_path.exists()

    treat.run(callback, max_gen=100)

    l5 = treat.replicates[5].get_lineage()
    del l5
    del treat

    def re_callback(rep, max_gen):
        lineage = rep.get_lineage()
        assert lineage.generation == max_gen

    treat = L.Treatment(path, p_3x2)
    treat.run(re_callback, max_gen=100)

