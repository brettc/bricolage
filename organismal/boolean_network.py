import logging
logging.basicConfig()
log = logging.getLogger("")


import itertools
import copy
import networkx as nx
import matplotlib.pyplot as plt

import numpy as np
randomizer = np.random


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


class Parameters(object):
    def __init__(self, **kwargs):
        # Defaults are provided here
        self.population_size = 20
        self.cis_count = 3
        self.cue_shapes = 2
        self.reg_shapes = 3
        self.out_shapes = 1
        self.ops = 'and or not_x not_y'.split()
        self.mutation_rate = .01

        self._override(kwargs)
        self._init()

    def _override(self, kwargs):
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                log.warning("'%s' is not a valid setting", k)

    def _init(self):
        self._calc_sizes()
        self._init_ops()
        self._init_envs()

    def random_sub(self):
        return randomizer.randint(*self.sub_range)

    def random_pub(self):
        # Only randomize the regulatory outputs
        return randomizer.randint(*self.reg_range)

    def random_op(self):
        i = randomizer.randint(0, self.ops_size)
        return self.ops_list[i]

    def _init_envs(self):
        self.environments = np.array([
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
        self.env_count = pow(2, self.cue_shapes)

        # This is per CIS module, rather than per site
        self.individual_mutation_rate = self.gene_count * self.cis_count * self.mutation_rate

        # Lastly, work out our char_size
        # Should never need more than 2, as 99 is really big
        if len(self.cue_signals) <= 10 and \
                len(self.reg_signals) <= 10 and \
                len(self.out_signals) <= 10:
            self.char_size = 2
        else:
            self.char_size = 3

    def name_for_signal(self, s):
        sz = self.char_size-1
        # We use "mathematical" number, starting at 1.0
        if s in self.cue_signals:
            return "E{0:0{1:}d}".format(s+1 - self.cue_range[0], sz)
        if s in self.reg_signals:
            return "T{0:0{1:}d}".format(s+1 - self.reg_range[0], sz)
        return "P{0:0{1:}d}".format(s+1 - self.out_range[0], sz)

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
        # TODO: Something sane here
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
        self.publish = p.pub_range[0] + i
        self.cis_modules = [CisModule(p) for _ in range(p.cis_count)]

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

# TODO: This needs a factory, to control the caching, mutation etc

class Network(object):
    def __init__(self, parent_or_params, generation=0, identifier=0):
        self.generation = generation
        self.identifier = identifier

        if isinstance(parent_or_params, Network):
            # This should only be called
            self._copy(parent_or_params)
            return
        
        params = parent_or_params
        self.params = params
        gc = self.params.gene_count
        self.genes = [Gene(self, i) for i in range(gc)]
        self.attractors = []
        self.rates = np.zeros((len(params.environments), params.out_shapes))
        self._init_states()

    def _copy(self, parent):
        # TODO: Move to factory -- 
        # Duplicate, but copy the Genes, rather than reference them
        self.params = parent.params
        self.genes = [copy.deepcopy(g) for g in parent.genes]

        # Don't update, this, as we'll be mutating
        self.attractors = []
        self.rates = np.zeros((len(self.params.environments), self.params.out_shapes))

    def mutated(self, number, generation=0, identifier=0):
        child = Network(self, generation, identifier)
        child.mutate(number)
        child._init_states()
        return child

    def _init_states(self):
        c = Cell(self.params)
        for i, e in enumerate(self.params.environments):
            c.set_environment(e)
            a = np.array(self.attractor(c))
            self.attractors.append(a)
            self.rates[i] = a.mean(axis=0)[-self.params.out_shapes:]

    def mutate(self, number):
        """We need gene / cis, then what to mutate? 
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
        self.env_size = params.cue_shapes
        self.products = params.product_array()
        self.activation = params.activation_array()

    def deactivate(self):
        self.activation[:] = 0

    def set_environment(self, env):
        self.products[:self.env_size] = env

    def decay(self):
        self.products[self.env_size:] = 0


class Target(object):
    def __init__(self, params, fn):
        self.params = params
        envs = self.params.environments
        self.opts = np.zeros((len(envs), params.out_shapes), dtype=float)
        for i, e in enumerate(envs):
            t = tuple(e)
            opt = fn(*t)
            self.opts[i] = opt


class Population(object):
    def __init__(self, params):
        self.params = params
        self.generation = 0
        self.networks = [Network(params, self.generation, i*self.generation) 
                         for i in range(params.population_size)]
        self.fitnesses = np.zeros(len(self.networks), dtype=float)

    def calc_fitness(self, target):
        for i, n in enumerate(self.networks):
            # Get the differences between desired and achieved
            diffs = abs(n.rates - target.opts)
            # Make them into fitness scores and normalise
            scores = (1 - diffs) / self.params.env_count

            # Sum across all challenges and scale by fitness contribution
            # summed = scores.sum(axis=1) * ch.fitness_contribution
            # summed = scores.sum(axis=1)

            fitness = scores.sum()
            fitness /= float(self.params.out_shapes)
            self.fitnesses[i] = fitness

    def roulette_selection(self):
        fit = self.fitnesses
        if np.any(fit < 0.0):
            log.error("Negative fitness in generation %s", self.generation)
            return 

        # Check there are some non-zero. Otherwise try again
        if np.all(fit == 0):
            log.warning("Generation %d: All fitnesses are zero, "
                        "so no selection", self.generation)
            # We do nothing, just let everything get mutated next time
            return

        # We calculate relative fitness by normalising the differences across
        # the population, with those that are alive (> 0.0)
        select_alive = np.where(fit > 0.0)
        alive = fit[select_alive]
        mn = min(alive)
        mx = max(alive)

        if mn == mx:
            log.warning("Generation %d: No variation in fitness "
                        "and no selection", self.generation)
            return

        # TODO: Maybe add a base amount (otherwise the very lowest get ejected
        relfit = (alive - mn) * 1.0 / (mx - mn)

        # We're pumping up the differences here. This makes it more likely that
        # small differences will be enough to get something into the next
        # generation. Should pbly make this into a parameter
        # TODO: Reinstate
        relfit = np.exp(relfit * 2.0) - 1.0

        fit[select_alive] = relfit

        # We'll fill them up by selecting from the parent population
        # Generate a list of increasing values to select from ...
        cumsum = fit.cumsum()
        # ... a list of random values to pick in the list
        sel_values = randomizer.uniform(0.0, cumsum[-1], len(self.networks))

        # Now return the corresponding indexes
        sel_indexes = np.searchsorted(cumsum, sel_values)

        # We keep these around things like attractors can do an update too
        self.selection_indexes = sel_indexes


    def new_generation(self):
        """Using the selection indexes, create a new generation"""
        new_networks = [self.networks[i] for i in self.selection_indexes]
        self.networks = new_networks

    def mutate(self):
        nmutations = randomizer.poisson(
            self.params.individual_mutation_rate, 
            self.params.population_size
        )

        for i in range(self.params.population_size):
            m = nmutations[i]
            if m:
                self.networks[i] = self.networks[i].mutated(m)


class SignalGraph(object):

    def signal_to_node(self, s):
        p = self.network.params
        if s in p.cue_signals:
            return 'E', s + 1 - p.cue_range[0]
        if s in p.reg_signals:
            return 'T', s+ 1 - p.reg_range[0]
        return 'P', s + 1 - p.out_range[0]

    def __init__(self, network):
        """Make a graph that shows signals connecting graphs"""
        self.network = network
        self.nx_graph = nx.DiGraph()

        for gene in self.network.genes:
            # node_name = gene.id_string
            print gene
            gene_node = 'G', gene.index

            out_node = self.signal_to_node(gene.publish)
            self.nx_graph.add_edge(gene_node, out_node)

            for c in gene.cis_modules:
                if c.active:
                    x_node = self.signal_to_node(c.x)
                    y_node = self.signal_to_node(c.y)
                    self.nx_graph.add_edge(x_node, gene_node)
                    self.nx_graph.add_edge(y_node, gene_node)

    def draw(self):
        # plt.ioff()
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)

        # Relabeling
        labels = {}
        for n in self.nx_graph.nodes_iter():
            if isinstance(n, tuple):
                labels[n] = n[0]+str(n[1])
            else:
                labels[n] = n
    
        # Draw the graph using graphvis layout
        pos = nx.graphviz_layout(self.nx_graph, prog='dot')
        nx.draw(self.nx_graph, pos=pos, ax=ax, node_size=500, font_size=6, node_color='w',
                with_labels=False)
        nx.draw_networkx_labels(self.nx_graph, pos, labels, font_size=8)
        # nx.draw_networkx_nodes(g,pos,nodelist=[largest_hub],node_size=300,node_color='r')
        #nx.draw_graphviz(g, ax=ax, fontsize=8, nodesize=100)

        # Stream it back as SVG
        output = open('x.png', 'wb')
        fig.savefig(output, format='png')
        # plt.clear() # This is the BUG I think
        # plt.ion() # turn interactive mode back on

if __name__ == '__main__':
    randomizer.seed(5)
    p = Parameters(
        cis_count=4,
        out_shapes=1,
        reg_shapes=3,
        cue_shapes=3,
        mutation_rate=.01,
        population_size=1000,
        ops='and or not_x not_y x_and_not_y y_and_not_x'.split(),
    )

    def ff(a, b, c):
        # return a and b, b and not a, a or b
        if a and b or c:
            return .5
        return 1.0

    t = Target(p, ff)
    pop = Population(p)
    n = pop.networks[0]
    g = SignalGraph(n)
    g.draw()
    for i in range(50):
        pop.calc_fitness(t)
        print max(pop.fitnesses)
        if max(pop.fitnesses) == 1.0:
            break
        pop.roulette_selection()
        pop.new_generation()
        pop.mutate()

    n = pop.networks[0]
    g = SignalGraph(n)
    g.draw()


    # n = Network(p)
    # c = n.mutated(5)
    
    # print n.describe()
    # print c.describe()
    # print n.attractors
    # print n.rates
    # c = Cell(p)
    # print n.describe()
    # c.expose(n)
    # n.mutate(2)
    # print '---------'
    # print n.describe()
    # c.expose(n)



