#!env python
import click

# This will initialise the folders.
from folders import data_dir
from bricolage import threshold3
from bricolage.lineage import SnapshotLineage


class GenerateError(Exception):
    pass


# A holder for the database.
_data = dict([
    (nm, (data_dir / nm).with_suffix('.db'))
    for nm in ['xor', 'bowtie']
])


# For external use
def get_database(name, readonly=False):
    dbpath = _data[name]
    assert dbpath.exists()
    return SnapshotLineage(dbpath, readonly=readonly)


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
        population_size=1000, mutation_rate=.002)

    with SnapshotLineage(dbpath, params=p) as lin:
        if not lin.targets:
            lin.add_target(bowtie_target)

        for g, b in select_till(lin, good_for=1000):
            click.secho("At generation {}, best is {}".format(g, b))


if __name__ == "__main__":
    generate()
