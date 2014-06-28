import logging
logging.basicConfig()
log = logging.getLogger("")

import numpy as np
randomizer = np.random
import itertools

def func_xor(x, y):
    return int((x or y) and not (x and y))

boolean_ops = {
    'not_x'       : ("~{0}"        , lambda x, y : not x)             ,
    'not_y'       : ("~{1}"        , lambda x, y : not y)             ,
    'and'         : ("{0} AND {1}" , lambda x, y : int(x and y))      ,
    'or'          : ("{0} OR {1}"  , lambda x, y : int(x or y))       ,
    'eq'          : ("{0} EQ {1}"  , lambda x, y : int(x == y))       ,
    'xor'         : ("{0} XOR {1}" , func_xor)                        ,
    'nand'        : ("{0} NAND {1}", lambda x, y : int(not (x and y))),
    'nor'         : ("{0} NOR {1}" , lambda x, y : int(not (x or y))) ,
    'x_and_not_y' : ("{0} AND ~{1}", lambda x, y : int(x and not y))  ,
    'y_and_not_x' : ("~{0} AND {1}", lambda x, y : int(y and not x))  ,
    'x_or_not_y'  : ("{0} OR ~{1}" , lambda x, y : int(x or not y))   ,
    'y_or_not_x'  : ("~{0} OR {1}" , lambda x, y : int(y or not x))   ,
}

class Defaults:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def apply(self, obj):
        obj.__dict__.update(self.__dict__)

_defaults = Defaults(
    # reg_gene_count=5,
    cis_count=3,
    cue_shapes=2,
    reg_shapes=3,
    out_shapes=1,
    ops='and or not_x not_y'.split(),
)

class Parameters(object):
    def __init__(self, **kwargs):
        _defaults.apply(self)
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                log.warning("'%s' is not a valid setting", k)
        self._init()

    def _init(self):
        self._calc_sizes()
        self._init_ops()

    def random_sub(self):
        return randomizer.randint(*self.sub_range)

    def random_pub(self):
        # Only randomize the regulatory outputs
        return randomizer.randint(*self.reg_range)

    def random_op(self):
        i = randomizer.randint(0, self.ops_size)
        return self.ops_list[i]

    def all_environments(self):
        return np.array([
            x for x in itertools.product(
                range(2), repeat=self.cue_shapes)])

    def _init_ops(self):
        self.ops_list = [boolean_ops[o] for o in self.ops]
        self.ops_size = len(self.ops_list)

    def product_array(self):
        return np.zeros(self.product_count, dtype=np.int8)

    def activation_array(self):
        return np.zeros(self.gene_count, dtype=np.int8)

    def _calc_sizes(self):
        # Calculate the total number of elements given the overlap
        # An example, with an overlap of 2, looks like this:
        # cue_range [ 0 1 2 3 ]
        # sub_range [ 0 1 2 3 4 5 ]
        # pub_range         [ 4 5 6 7 8 9 10 ]
        # out_range             [ 6 7 8 9 10 ]
        #
        # reg_range         [ 4 5 ] 
        self.reg_gene_count = self.reg_shapes

        self.sub_shapes = self.cue_shapes + self.reg_shapes

        self.product_count = self.cue_shapes + self.reg_shapes + \
            self.out_shapes

        self.sub_range = 0, self.sub_shapes
        self.reg_range = self.cue_shapes, self.sub_shapes
        self.pub_range = self.cue_shapes, self.sub_shapes + self.out_shapes
        self.cue_range = 0, self.cue_shapes
        self.out_range = self.sub_shapes, self.product_count

        # Let's make these, they'll be useful
        self.cue_signals = range(*self.cue_range)
        self.reg_signals = range(*self.reg_range)
        self.out_signals = range(*self.out_range)

        self.sub_signals = range(*self.sub_range)
        self.pub_signals = range(*self.pub_range)

        # Total gene count requires 
        self.gene_count = self.reg_gene_count + self.out_shapes


class CisModule(object):
    def __init__(self, params):
        self.active = randomizer.uniform(0, 1) > .5
        self.x = params.random_sub()
        self.y = params.random_sub()
        self.description, self.operation = params.random_op()

    def bind(self, cell):
        if self.active:
            a = cell.products[self.x]
            b = cell.products[self.y]
            if self.operation(a, b):
                return True
        return False

    def mutate(self, params):
        # Decide what to do...
        q = randomizer.uniform(0, 100)
        if q > 80:
            self.description, self.operation = params.random_op()
        elif q > 40:
            self.x = params.random_sub()
        else:
            self.y = params.random_sub()

        # Flip
        if randomizer.uniform(0, 1) > .5:
            self.active = not self.active

    def describe(self):
        return self.description.format(self.x, self.y)


class Gene(object):
    def __init__(self, network, i):
        self.index = i
        p = network.params
        self.cis_modules = [CisModule(p) for _ in range(p.cis_count)]
        
        # out_index = i - p.reg_gene_count
        # if out_index < 0:
        #     self.publish = p.random_pub()
        # else:
            # self.publish = p.out_signals[out_index]
        self.publish = p.pub_range[0] + i

    def bind(self, cell):
        for m in self.cis_modules:
            if m.bind(cell):
                cell.activation[self.index] = 1
                # One is sufficient
                return

    def transcribe(self, cell):
        if cell.activation[self.index]:
            cell.products[self.publish] = 1

    def describe(self):
        ds = [m.describe() for m in self.cis_modules if m.active]
        return '(' + ') OR ('.join(ds) + ') => {}'.format(self.publish)


class Network(object):
    def __init__(self, params):
        self.params = params
        gc = self.params.gene_count
        self.genes = [Gene(self, i) for i in range(gc)]


    def mutate(self, number):
        """
        We need gene / cis, then what to mutate? 
        """
        positions = self.params.gene_count * self.params.cis_count
        mutations = randomizer.randint(0, positions, number)
        for m in mutations:
            gene_pos, cis_pos = divmod(m, self.params.cis_count)
            cis = self.genes[gene_pos].cis_modules[cis_pos]
            cis.mutate(self.params)

    def cycle(self, cell):
        """Run a cycle through the network, updating the elements
        
        This is the core update. Everything counts on this working. We don't
        go for speed here.  The speedy version is in cython. But we make sure
        that this method and the cython method give the same results.
        """

        # This binds and activates
        cell.deactivate()
        for g in self.genes:
            g.bind(cell)

        # "Decay" everything not from the environment
        cell.decay()

        for g in self.genes:
            g.transcribe(cell)

    def attractor(self, cell):
        # Put starter position into transient
        seen = []

        while True:
            self.cycle(cell)
            cell_state = cell.products.copy()
            for i, s in enumerate(seen):
                if np.all(s == cell_state):
                    return seen[i:]
                    # transient = seen[:i]
            seen.append(cell_state)

        raise RuntimeError

    def expression_rate(self, cell):
        a = np.array(self.attractor(cell))
        return a.mean(axis=0)[-self.params.out_shapes:]

    def describe(self):
        gs = [g.describe() for g in self.genes]
        print '\n'.join(gs)
            

class Cell(object):
    def __init__(self, params):
        self.params = params
        self.env_size = self.params.cue_shapes
        self.products = params.product_array()
        self.activation = params.activation_array()

    def deactivate(self):
        self.activation[:] = 0

    def set_environment(self, env):
        self.env = env
        self.products[:self.env_size] = env

    def decay(self):
        self.products[self.env_size:] = 0

    def expose(self, n):
        for e in self.params.all_environments():
            c.set_environment(e)
            print e, n.expression_rate(c)
            # print '--'
            # # a = n.attractor(c)
            # # for s in a:
            # #     print s
            # print '-----'

if __name__ == '__main__':
    p = Parameters(
        cis_count=4,
        out_shapes=3,
        reg_shapes=5,
        cue_shapes=2,
        ops='and not_x not_y x_and_not_y y_and_not_x'.split(),
    )
    n = Network(p)
    c = Cell(p)
    print n.describe()
    c.expose(n)
    n.mutate(2)
    print '---------'
    print n.describe()
    c.expose(n)



