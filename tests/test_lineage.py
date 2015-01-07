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

# def test_numpy_export(c_3x2):
#     p1 = T.Population(c_3x2, 1000)
#     print p1
#     as_array = c_3x2.to_numpy(p1, mutations_only=True)
#     print as_array
#     p1.mutate(.0001)
#     print len(p1.mutated)
#     as_array = c_3x2.to_numpy(p1, mutations_only=True)
#     print len(as_array)
#     print as_array
#

    # assert as_array.dtype == c_3x2.dtype()
    # print c_3x2.dtype()
    #
    # p2 = T.Population(c_3x2, 0)
    # c_3x2.from_numpy(as_array, p2)
    #

def assert_pops_equal(p1, p2):
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
    fname = pathlib.Path(str(tmpdir)) / 'creation.db'
    params = T.Parameters(population_size=10)
    a = L.Lineage(fname, params, T.Constructor)
    del a
    b = L.Lineage(fname)
    del b

def test_restarting(tmpdir, p_3x2, target_3x2):
    fname = pathlib.Path(str(tmpdir)) / 'selection.db'

    # Set it all up
    a = L.Lineage(fname, p_3x2, T.Constructor)
    sel = T.SelectionModel(a.world)
    target = T.Target(a.world, target_3x2)

    # Run selection for 100 generations
    for i in range(100):
        a.next_generation(.01, target, sel)
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
        c.next_generation(.01, target, sel)
    del c

    # Relead yet again. Reload the 100th generation again. Check it is exactly
    # the same
    d = L.Lineage(fname)
    p2 = d.get_generation(100)
    z = [net.identifier for net in p2]
    assert y == z

    # This means everything is the same!
    assert_pops_equal(p1, p2)
    

