import logging
logging.basicConfig()
log = logging.getLogger("")


import itertools
import StringIO
import networkx as nx
import matplotlib.pyplot as plt

import numpy as np
randomizer = np.random
import pandas as pd
int_type = np.int
tiny_type = np.int8

class Parameters(object):
    def __init__(self, **kwargs):
        # Defaults are provided here
        self.reg_gene_count = 2
        self.reg_shapes = 2

        self.cue_shapes = 1
        self.out_shapes = 1
        self.ps1 = .5

        self.population_size = 1000
        self.mutation_rate = .01
        self.generations = 100

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
        self._init_envs()
        self._calc_mutator()

    def _init_envs(self):
        self.environments = np.array([
            x for x in itertools.product(
                range(2), repeat=self.cue_shapes)])

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
        #

        self.struct_gene_count = self.out_shapes
        self.gene_count = self.reg_gene_count + self.out_shapes

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

        self.reg_set = set(self.reg_signals)
        self.sub_set = set(self.sub_signals)

        # Total gene count requires 
        self.gene_count = self.reg_gene_count + self.out_shapes
        self.env_count = pow(2, self.cue_shapes)

        # Lastly, work out our char_size
        # Should never need more than 2, as 99 is really big
        if len(self.cue_signals) <= 10 and \
                len(self.reg_signals) <= 10 and \
                len(self.out_signals) <= 10:
            self.char_size = 2
        else:
            self.char_size = 3

        self.reg_gene_size = 3
        self.struct_gene_size = 2

        self.dna_size = self.reg_gene_count * self.reg_gene_size \
            + self.struct_gene_count * self.struct_gene_size

        self.fitness_contribution = np.array([[self.ps1], [1.0-self.ps1]])

    def _calc_mutator(self):
        # A rough but simple way to generate mutations in the dna. We just
        # select from this list, it gives us all the info
        mutate_info = []
        codon = 0
        strat_size = len(Gene.strategies)
        for i in range(self.reg_gene_count):
            mutate_info.append((codon, set(range(strat_size))))
            codon +=1 
            mutate_info.append((codon, self.sub_set))
            codon +=1 
            mutate_info.append((codon, self.reg_set))
            codon +=1 
        for i in range(self.struct_gene_count):
            mutate_info.append((codon, set(range(strat_size))))
            codon +=1 
            mutate_info.append((codon, self.reg_set))
            codon +=1 

        self.mutate_info = mutate_info
        self.mutate_size = len(mutate_info)

    def name_for_signal(self, s):
        sz = self.char_size-1
        # We use "mathematical" number, starting at 1.0
        if s in self.cue_signals:
            return "E{0:0{1:}d}".format(s+1 - self.cue_range[0], sz)
        if s in self.reg_signals:
            return "T{0:0{1:}d}".format(s+1 - self.reg_range[0], sz)
        return "P{0:0{1:}d}".format(s+1 - self.out_range[0], sz)


class Gene(object):
    # strategies = [(0, 1), (1, 0)]
    strategies = [(0, 0), (0, 1), (1, 0), (1, 1)]

    def __init__(self, network, i, dna):
        self.index = i
        self.dna = dna
        self.strat = dna[0]
        self.on, self.off = self.strategies[self.strat]
        self.subscribe = dna[1]
        if len(dna) == network.params.reg_gene_size:
            self.publish = dna[2]
        else:
            assert len(dna) == network.params.struct_gene_size
            p = network.params
            # Publish the 
            self.publish = i - p.reg_gene_count + p.out_range[0]

    def bind(self, cell):
        state = cell.products[self.subscribe]
        if state:
            cell.activation[self.index] = self.on
        else:
            cell.activation[self.index] = self.off

    def transcribe(self, cell):
        if cell.activation[self.index]:
            cell.products[self.publish] = 1

    def describe(self):
        return "Gene<s:{0.subscribe}, p:{0.publish}, d:{0.off}/{0.on}>".format(self)

    def __repr__(self):
        return self.describe()

    def __str__(self):
        return self.describe()


# TODO: This needs a factory, to control the caching, mutation etc
class Network(object):
    def __init__(self, params, dna):
        self.params = params
        self.genes = [Gene(self, i, d) for i, d in enumerate(dna)]

        # Flatten the dna
        self.dna = np.zeros(params.dna_size, dtype=int)
        i = 0
        for gene in dna:
            for site in gene:
                self.dna[i] = site
                i += 1

        self.attractors = []
        self.rates = np.zeros((len(params.environments), params.out_shapes))
        self._init_states()

        self.neighbours = set()

    def is_neighbour(self, other):
        diffs = sum(self.dna != other.dna)
        if diffs == 1:
            return True
        return False

    def _init_states(self):
        c = Cell(self.params)
        for i, e in enumerate(self.params.environments):
            c.set_environment(e)
            a = np.array(self.attractor(c))
            self.attractors.append(a)
            self.rates[i] = a.mean(axis=0)[-self.params.out_shapes:]

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


class AllPossible(object):
    def __init__(self, params, target):
        self.params = params
        self.target = target
        self.networks = []
        self.fitnesses = []
        self.pathlen = []
        self._generate()
        self._calc_fitness()

    def _generate(self):
        reg_genes = [self._gene_poss()] * (self.params.reg_gene_count)
        struct_genes = [self._gene_poss(reg=False)] * (
            self.params.gene_count - self.params.reg_gene_count)
        genes = reg_genes + struct_genes
        for dna in itertools.product(*genes):
            n = Network(self.params, dna)
            self.networks.append(n)

    def _gene_poss(self, reg=True):
        all_pub = self.params.reg_signals
        all_sub = self.params.sub_signals
        on_off = range(len(Gene.strategies))
        if reg:
            generate_params = [on_off, all_sub, all_pub]
        else:
            generate_params = [on_off, all_pub]
        return [_ for _ in itertools.product(*generate_params)]

    def _calc_fitness(self):
        t = self.target
        self.fitnesses = []
        for i, n in enumerate(self.networks):
            # Get the differences between desired and achieved
            diffs = abs(n.rates - t.opts)
            # Make them into fitness scores and normalise
            # scores = (1 - diffs) / self.params.env_count
            scores = (1 - diffs) * self.params.fitness_contribution

            # print diffs
            # print scores
            # print scores2
            # print '--'
            fitness = scores.sum()

            # fitness = scores.sum()
            fitness /= float(self.params.out_shapes)
            self.fitnesses.append(fitness)

            if fitness == 1.0:
                l = self._calc_path(n)
            else:
                l = None
            self.pathlen.append(l)

            # Shitty, but this will be useful
            n.pathlen = l
            n.fitness = fitness

    def _calc_path(self, n):
        n.graph = SignalGraph(n)
        plen = len(nx.shortest_path(n.graph.nx_graph, ('C',1), ('A',1)))
        plen = (plen-1)/2
        return plen

    def get_dataframe(self):
        df = pd.DataFrame(dict(
            network=self.networks, 
            fitness=self.fitnesses,
            pathlen=self.pathlen,
            ))

        _add_strategy_cols(df)
        return df
        
def _add_strategy_cols(df):
    nets = df['network']
    cols = {}
    for n in nets:
        for i, g in enumerate(n.genes):
            cols.setdefault('g%d-Cpub' % i, []).append(g.publish)
            cols.setdefault('g%d-Bsub' % i, []).append(g.subscribe)
            cols.setdefault('g%d-Astrat' % i, []).append(g.dna[0])

    keys = cols.keys()
    keys.sort()
    for k in keys:
        df[k] = cols[k]

class Neighbourhood(object):
    def __init__(self, networks):
        self.G = nx.Graph()
        self.networks = list(networks)
        self._calc_neighbours()
        self._calc_subgraphs()

    def _calc_neighbours(self):
        networks = self.networks
        ncount = len(networks)
        for i, this_net in enumerate(self.networks):
            self.G.add_node(this_net)
            for j in range(i+1, ncount):
                other_net = networks[j]
                if this_net.is_neighbour(other_net):
                    self.G.add_edge(this_net, other_net)

                    this_net.neighbours.add(other_net)
                    other_net.neighbours.add(this_net)

    def _calc_subgraphs(self):
        self.components = nx.connected_components(self.G)
        for i, comp in enumerate(self.components):
            for net in comp:
                net.subgraph = i

    def get_dataframe(self):
        nets = self.networks
        df = pd.DataFrame(dict(
            network=self.networks, 
            fitness=[n.fitness for n in nets],
            pathlen=[n.pathlen for n in nets],
            component=[n.subgraph for n in nets],
            edges=[len(self.G.edges(n)) for n in nets],
            ))

        _add_strategy_cols(df)
        return df


class SignalGraph(object):

    def signal_to_node(self, s):
        p = self.network.params
        if s in p.cue_signals:
            return 'C', s + 1 - p.cue_range[0]
        if s in p.reg_signals:
            return 'S', s + 1 - p.reg_range[0]
        return 'A', s + 1 - p.out_range[0]

    def __init__(self, network):
        """Make a graph that shows signals connecting graphs"""
        self.network = network
        self.nx_graph = nx.DiGraph()

        for gene in self.network.genes:
            # gene_node = 'G', gene.index+1
            gene_node = gene

            out_node = self.signal_to_node(gene.publish)
            self.nx_graph.add_edge(gene_node, out_node)

            in_node = self.signal_to_node(gene.subscribe)
            self.nx_graph.add_edge(in_node, gene_node)

    # def draw(self):
    def _repr_svg_(self):
        # plt.ioff()
        fig = plt.figure(figsize=(15,4))
        ax = fig.add_subplot(111)

        # Relabeling
        labels = {}
        for n in self.nx_graph.nodes_iter():
            if isinstance(n, tuple):
                labels[n] = n[0]+str(n[1])
            else:
                labels[n] = "G{}-{}{}".format(n.index+1, n.off, n.on)
    
        # Draw the graph using graphvis layout
        # pos = nx.fruchterman_reingold_layout(self.nx_graph)
        pos = nx.graphviz_layout(self.nx_graph, prog='dot', args='-Grankdir=LR')
        nx.draw(self.nx_graph, pos=pos, ax=ax, node_size=800, font_size=6, node_color='w',
                with_labels=False)
        nx.draw_networkx_labels(self.nx_graph, pos, labels, font_size=8)
        # nx.draw_networkx_nodes(g,pos,nodelist=[largest_hub],node_size=300,node_color='r')
        #nx.draw_graphviz(g, ax=ax, fontsize=8, nodesize=100)

        # Stream it back as SVG
        output = StringIO.StringIO()
        fig.savefig(output, format='svg')
        # plt.clear() # This is the BUG I think
        plt.ion() # turn interactive mode back on
        # return output.getvalue()
        #
        
class Simulation(object):
    def __init__(self, params, allposs):
        self.params = params
        self.allposs = allposs
        self.networks_by_dna = {}
        for n in allposs.networks:
            self.networks_by_dna[tuple(n.dna)] = n
        self.population = []

    def run(self, seedn):
        self.population = [seedn] * self.params.population_size
        for g in range(self.params.generations):
            self.mutate()
            self.select()

    def select(self):
        fit = np.array([n.fitness for n in self.population])
        cumsum = fit.cumsum()
        # ... a list of random values to pick in the list
        sel_values = randomizer.uniform(0.0, cumsum[-1], self.params.population_size)

        # Now return the corresponding indexes
        sel_indexes = np.searchsorted(cumsum, sel_values)

        new_pop = [self.population[i] for i in sel_indexes]
        self.population = new_pop

    def mutate(self):
        mrate = self.params.mutation_rate
        nmutations = randomizer.poisson(
            mrate * self.params.population_size * self.params.mutate_size)
        indivs = randomizer.randint(0, self.params.population_size, nmutations)
        for i in indivs:
            n = self.population[i]
            new_n = self.mutated(n)
            self.population[i] = new_n

    def mutated(self, n):
        mi = self.params.mutate_info
        assert len(mi) == len(n.dna)
        pos = randomizer.randint(0, len(mi))
        codon, poss_vals = mi[pos]

        # Ignore the current value
        current_val = n.dna[codon]
        new_vals = poss_vals.copy()
        new_vals.remove(current_val)
        new_vals = list(new_vals)

        # Might be only one other choice
        if len(new_vals) == 1:
            new_val = new_vals[0]
        else:
            idx = randomizer.randint(len(new_vals))
            new_val = new_vals[idx]

        new_dna = n.dna.copy()
        new_dna[codon] = new_val

        return self.networks_by_dna[tuple(new_dna)]


if __name__ == '__main__':
    # randomizer.seed(5)
    p = Parameters(
        reg_gene_count=2,
        out_shapes=1,
        reg_shapes=2,
        cue_shapes=1,
        ps1=.1,
        generations=200,
    )

    def ff(a):
        # return a and b, b and not a, a or b
        return a

    t = Target(p, ff)
    x = AllPossible(p, t)
    s = Simulation(p, x)
    s.run(x.networks[0])
    for n in s.population:
        print n.pathlen, n.fitness,

    # n = pop.networks[0]
    # g = SignalGraph(n)
    # g.draw()
    # for i in range(20):
    #     pop.calc_fitness(t)
    #     print pop.fitnesses.mean()
    #     pop.roulette_selection()
    #     pop.new_generation()
    #     pop.mutate()
    #
    #
