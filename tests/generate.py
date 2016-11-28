#!env python
import click
import numpy
import itertools

from folders import data_dir
from bricolage import threshold3, logic2
from bricolage.lineage import SnapshotLineage
from bricolage.core import InputType, ScoringMethod
from bricolage.core_ext import SelectionModel
from bricolage.targets_ext import NoisyTarget, DefaultTarget


class GenerateError(Exception):
    pass


# A holder for the database.
_data = dict([
    (nm, (data_dir / nm).with_suffix('.db'))
    for nm in ['xor', 'bowtie', 'perturb', 'three', 'double_bow']
])

_numpy_dumps = dict([
    (nm, (data_dir / nm).with_suffix('.npy'))
    for nm in ['noisy']
])


# For external use
def get_database(name, readonly=False):
    dbpath = _data[name]
    assert dbpath.exists()
    return SnapshotLineage(dbpath, readonly=readonly)

def get_numpy_dump(name, readonly=False):
    pth = _numpy_dumps[name]
    assert pth.exists()
    return numpy.load(str(pth))

@click.group(chain=True)
def generate():
    pass


def select_till(lin, good_for=1, return_every=100):
    got_1 = 0
    while 1:
        lin.next_generation()
        w, b = lin.population.worst_and_best()

        if lin.generation % 100 == 0:
            yield lin.generation, b

        if b == 1.0:
            got_1 += 1

        # Run on for 1000 extra generations
        if got_1 == good_for:
            break


def get_dbpath(name, overwrite):
    global _data
    dbpath = _data[name]

    if dbpath.exists() and not overwrite:
        click.secho("{} already exists".format(str(dbpath)))
        return None

    return dbpath


def xor_target(a, b):
    if (a or b) and not (a and b):
        return [.5, 1]
    return [0, 0]


@generate.command()
@click.option('--overwrite', is_flag=True, default=False)
def xor(overwrite):
    dbpath = get_dbpath('xor', overwrite)
    if dbpath is None:
        return

    p = threshold3.Parameters(
        seed=3, cis_count=2, reg_channels=4, out_channels=2, cue_channels=2,
        population_size=1000, mutation_rate=.001,)

    with SnapshotLineage(dbpath, params=p) as lin:
        if not lin.targets:
            lin.add_target(xor_target)

        for g, b in select_till(lin, good_for=1000):
            click.secho("At generation {}, best is {}".format(g, b))


def bowtie_target(a, b, c):
    if (a and not c) or (b and c):
        return [1, .5, .25]
    return [0, 0, 0]


@generate.command()
@click.option('--overwrite', is_flag=True, default=False)
def bowtie(overwrite):
    dbpath = get_dbpath('bowtie', overwrite)
    if dbpath is None:
        return

    p = threshold3.Parameters(
        seed=8, cis_count=2, reg_channels=8, out_channels=3, cue_channels=3,
        population_size=1000, mutation_rate=.002,
        input_type = InputType.PULSE,
    )

    with SnapshotLineage(dbpath, params=p) as lin:
        if not lin.targets:
            lin.add_target(bowtie_target)

        for g, b in select_till(lin, good_for=1000):
            click.secho("At generation {}, best is {}".format(g, b))

def perturb_target(a, b):
    if (a or b) and not (a and b):
        return 1
    return 0


@generate.command()
@click.option('--overwrite', is_flag=True, default=False)
def perturb(overwrite):
    dbpath = get_dbpath('perturb', overwrite)
    if dbpath is None:
        return

    p = threshold3.Parameters(
        seed=1, cis_count=2, reg_channels=4, out_channels=1, cue_channels=2,
        population_size=1000, mutation_rate=.002, input_type=InputType.PULSE,
        target_class=NoisyTarget,
    )

    with SnapshotLineage(dbpath, params=p) as lin:
        if not lin.targets:
            lin.targets.append(NoisyTarget(lin.world, perturb_target,
                                           perturb_count=3, perturb_prop=.2))
            # lin.targets.append(DefaultTarget(lin.world, perturb_target))
            lin.set_target(0)

        for g, b in select_till(lin, good_for=10):
            click.secho("At generation {}, best is {}".format(g, b))


def make_mapping():
    cats = [0, 0, 0, 1, 1, 1, 2, 2]
    assert len(cats) == 8
    perms = itertools.permutations(cats)

    # Mix it up
    for i in range(5):
        chosen = perms.next()

    # All possible inputs
    inputs = [_ for _ in itertools.product([0, 1], repeat=3)]
    mapping = dict([(a, b) for a, b in zip(inputs, chosen)])
    return mapping


def three_target(a, b, c):
    mapping = make_mapping()
    cat = mapping[(a, b, c)]
    if cat == 0:
        return 0, 0, 1, 1
    if cat == 1:
        return 1, 1, 0, 0
    return 0, 1, 1, 0


@generate.command()
@click.option('--overwrite', is_flag=True, default=False)
def three(overwrite):
    dbpath = get_dbpath('three', overwrite)
    if dbpath is None:
        return

    # Good values seed=2
    p = threshold3.Parameters(
        seed=15, cis_count=3, reg_channels=8, out_channels=4, cue_channels=3,
        # seed=2, cis_count=3, reg_channels=8, out_channels=4, cue_channels=3,
        population_size=1000, mutation_rate=.002, input_type=InputType.PULSE,
    )

    with SnapshotLineage(dbpath, params=p) as lin:
        if not lin.targets:
            lin.targets.append(NoisyTarget(lin.world, three_target,
                                           perturb_count=1, perturb_prop=.1))
            # lin.targets.append(DefaultTarget(lin.world, three_target))
            lin.set_target(0)

        for g, b in select_till(lin, good_for=1000):
            click.secho("At generation {}, best is {}".format(g, b))


def xor_func(a, b):
    return (a or b) and not (a and b)

def fit_func1(a, b):
    return xor_func(a, b)


def make_noisy():
    p = threshold3.Parameters(
        seed=1, cis_count=4, reg_channels=8, out_channels=3, cue_channels=3,
    )
    world = threshold3.World(p)
    factory = threshold3.Factory(world)
    ntarget = NoisyTarget(world, bowtie_target, perturb_count=3, perturb_prop=.2)
    pop = threshold3.Population(factory, 1000)
    select = SelectionModel(world)
    for i in range(20):
        pop.select(select)
        pop.mutate(.002)
        pop.assess(ntarget)

    fits = pop.fitnesses
    return fits


@generate.command()
@click.option('--overwrite', is_flag=True, default=False)
def noisy_fitnesses(overwrite):
    numpy.save(str(_numpy_dumps['noisy']), make_noisy())

def double_bow_target(a, b, c, d, e, f):
    ret = []
    if (a and not b) or (b and not c):
        ret.extend([0, 0, 1, 1])
    else:
        ret.extend([1, 1, 0, 0])

    if (d and not e) or (d and not f) or (not e and not f):
        ret.extend([0, 0, 1, 1])
    else:
        ret.extend([1, 1, 0, 0])
    return ret


@generate.command()
@click.option('--overwrite', is_flag=True, default=False)
def double_bow(overwrite):
    dbpath = get_dbpath('double_bow', overwrite)
    if dbpath is None:
        return

    p = logic2.Parameters(
        seed=21166, cis_count=3, reg_channels=14, out_channels=8, cue_channels=6,
        population_size=5000, mutation_rate=.001,
    )

    with SnapshotLineage(dbpath, params=p) as lin:
        lin.add_target(DefaultTarget(lin.world, double_bow_target,
                       scoring_method=ScoringMethod.EXPONENTIAL_VEC,
                       strength=.2))

        for g, b in select_till(lin, good_for=1000):
            click.secho("At generation {}, best is {}".format(g, b))

if __name__ == "__main__":
    generate()
