#!env python
"""
Generate data for the tests. 

Usage:
  generate all [--overwrite] 
  generate bowtie [--overwrite] 
  generate xor [--overwrite] 
  generate -h | --help 

Options:
  --overwrite       Overwrite the current data if it exists
  -h --help         Show this screen.
"""
__version__ = '0.0.1'

class GenerateError(Exception):
    pass


from docopt import docopt
import sys
from clint.textui import puts, indent, colored

# This initialise folders 
from folders import data_dir
from bricolage import threshold3
from bricolage.lineage import SnapshotLineage

_commands = {}
_data = {}

# For external use
def get_database(name, readonly=False):
    dbpath = _data[name]
    assert dbpath.exists()
    return SnapshotLineage(dbpath, readonly=readonly)

def command(func):
    nm = func.__name__
    dbpath = (data_dir / nm).with_suffix('.db')
    strpath = str(dbpath)

    def wrapped(arguments):
        puts(colored.blue("Beginning command {} ...".format(nm)))
        with indent(4):
            if dbpath.exists() and not arguments.options.overwrite:
                puts(colored.red("Database {} already exists".format(strpath)))
                return
            puts(colored.green("Using database {}".format(strpath)))
            func(dbpath, arguments)
        puts(colored.blue("End command {} ...".format(nm)))

    _commands[nm] = wrapped
    _data[nm] = dbpath
    return wrapped


def select_till(L, good_for=1, return_every=100):
    got_1 = 0
    while 1:
        L.next_generation()
        w, b = L.population.worst_and_best()

        if L.generation % 100 == 0:
            yield L.generation, b

        if b == 1.0:
            got_1 += 1

        # Run on for 1000 extra generations
        if got_1 == good_for:
            break

def xor_target(a, b):
    if (a or b) and not (a and b):
        return [.5, 1]
    return [0, 0]

@command
def xor(dbpath, arguments):
    p = threshold3.Parameters(
        seed=3, cis_count=2, reg_channels=4, out_channels=2, cue_channels=2,
        population_size=1000, mutation_rate=.001,)

    with SnapshotLineage(dbpath, params=p) as L:
        if not L.targets:
            L.add_target(xor_target)

        for g, b in select_till(L, good_for=1000):
            puts("At generation {}, best is {}".format(g, b))

def bowtie_target(a, b, c):
    if (a and not c) or (b and c):
        return [1, .5, .25]
    return [0, 0, 0]

@command
def bowtie(dbpath, arguments):
    p = threshold3.Parameters(
        seed=8, cis_count=2, reg_channels=8, out_channels=3, cue_channels=3,
        population_size=1000, mutation_rate=.002)
    
    with SnapshotLineage(dbpath, params=p) as L:
        if not L.targets:
            L.add_target(bowtie_target)

        for g, b in select_till(L, good_for=1000):
            puts("At generation {}, best is {}".format(g, b))

def main(arguments):
    if not data_dir.exists():
        puts("Making data folder {}".format(str(data_dir)))
        data_dir.mkdir()

    try:
        for cmd in arguments.current_commands:
            func = _commands[cmd]
            func(arguments)
    except:
        raise
    return 0

class Arguments(object):
    """Wrapper for docopt arguments -- makes processing much simpler"""
    class Container(object):
        """Empty class allows referencing without dictionary fluff"""
        pass

    def __init__(self, arguments):
        all_cmds = []
        current = None
        self.options = Arguments.Container()
        self.arguments = Arguments.Container()

        # Turn arguments of different kinds into sensible stuff
        for arg, value in arguments.items():
            if arg.startswith('--'):
                setattr(self.options, arg[2:], value)
            elif arg.startswith('<'):
                # TODO: Move out sanitisation
                new_arg = arg.replace("<", "").replace(">", "")
                setattr(self.arguments, new_arg, value)
            else:
                if arg != 'all':
                    all_cmds.append(arg)
                if arguments[arg]:
                    current = arg
        self.all_commands = all_cmds
        if current == 'all':
            self.current_commands = all_cmds
        else:
            self.current_commands = [current]

if __name__ == "__main__":
    arguments = Arguments(docopt(__doc__, version=__version__))
    sys.exit(main(arguments))
