"""High level analysis tools
"""
import logtools
import numpy as np
from .analysis_ext import (MutualInfoAnalyzer, OutputControlAnalyzer, RelevantControlAnalyzer)
from .analysis import AverageControlAnalyzer
from .lineage import FullLineage
from .experimentdb import StatsGroupRecord, StatsRecord, StatsReplicateRecord
from .neighbourhood import PopulationNeighbourhood
from bricolage.experiment import Experiment

log = logtools.get_logger()


class StatsVisitor(object):
    def __init__(self, experiment, stats_classes):
        assert isinstance(experiment, Experiment)
        self.experiment = experiment
        self.session = experiment.database.session
        self.replicate = None
        self.calc_count = 0

        self.todo = set()
        # Make sure there are not repeats
        for kls in set(stats_classes):
            k = kls()
            # We need a tag attribute.
            assert hasattr(k, 'tag')
            self.todo.add(k)

        # Load all stats groups already done.
        # NOTE: Not sure this is the best way to do it yet...
        self.done = {}

        log.info("Loading all stats groups...")
        for grp in self.session.query(StatsGroupRecord).all():
            self.done[(grp.treatment_id, grp.replicate_id, grp.generation, grp.tag)] = grp
        log.info("... finished loading.")

    def is_done(self, rep, gen_num, tag):
        return (rep.treatment.seq, rep.seq, gen_num, tag) in self.done

    def add_stats_group(self, rep, gen_num, tag):
        self.session.add(StatsGroupRecord(rep, gen_num, tag))

    def visit_lineage(self, rep, lin):
        self.replicate = rep
        for stats in self.todo:
            stats.init_lineage(rep, lin)

    def wants_generation(self, gen_num):
        for stats in self.todo:
            if not self.is_done(self.replicate, gen_num, stats.tag):
                # Ok -- at least one thing is missing
                return True
        return False

    def visit_generation(self, gen_num, pop):
        log.info("Doing stats for generation {}, {}".format(self.replicate.path, gen_num))
        for stats in self.todo:
            if self.is_done(self.replicate, gen_num, stats.tag):
                log.info("Skipping already created group for %s, %s", gen_num, stats.tag)
                continue

            self.add_stats_group(self.replicate, gen_num, stats.tag)
            named_values = stats.calc_stats(pop)
            self.session.add_all([
                StatsRecord(
                    self.replicate,
                    gen_num,
                    "{}_{}".format(stats.tag, n), v, stats.tag)
                for n, v in named_values
            ])

        self.session.commit()

    #     if self.calc_count > 20:
    #         log.info("Committing to database...")
    #         self.session.commit()
    #         self.calc_count = 0
    #     else:
    #         self.calc_count += 1
    #
    # def leave_lineage(self, gen_num, pop):
    #     log.info("Committing to database...")
    #     self.session.commit()


class StatsAverageControl(object):
    tag = "AC"

    def __init__(self):
        self.analyzer = None
        self.regs = None

    def init_lineage(self, rep, lin):
        self.analyzer = AverageControlAnalyzer(lin.world)
        self.regs = lin.params.reg_channels
        self.out = lin.params.out_channels

    def calc_stats(self, pop):
        ai = self.analyzer.calc_info(pop)

        # Summarize across the population
        control = ai.control.mean(axis=0)
        entropy = ai.entropy.mean(axis=0)

        vals = []

        # Record the mean of all information measures
        regs = self.regs
        out = self.out
        for i, c in enumerate(range(regs)):
            for j in range(out):
                vals.append(('C_{}_{}'.format(c + 1, j + 1), control[i, j]))
                vals.append(('E_{}_{}'.format(c + 1, j + 1), entropy[i, j]))

        vals.extend([
            ('C_MEAN', control.mean()),
            ('E_MEAN', entropy.mean()),
        ])
        return vals


class StatsOutputControl(object):
    tag = "OC"

    def __init__(self):
        self.analyzer = None

    def init_lineage(self, rep, lin):
        targ = lin.targets[0]
        tset = targ.calc_distinct_outputs()
        self.analyzer = OutputControlAnalyzer(lin.world, tset)
        self.regs = lin.params.reg_channels

    def calc_stats(self, pop):
        ai = np.asarray(self.analyzer.analyse_collection(pop))
        assert not np.any(np.isnan(ai.ravel()))

        # Summarize across the population
        ameans = np.mean(ai, axis=0)

        vals = []

        # Record the mean of all information measures
        regs = self.regs
        for i, c in enumerate(range(regs)):
            vals.append(('C_{}'.format(c + 1), ameans[i, 0]))
            vals.append(('E_{}'.format(c + 1), ameans[i, 1]))
            vals.append(('W_{}'.format(c + 1), ameans[i, 2]))

        reg_mean = ameans.mean(axis=0)
        vals.extend([
            ('C_MEAN', reg_mean[0]),
            ('E_MEAN', reg_mean[1]),
            ('W_MEAN', reg_mean[2]),
            ('C_MAX', ai[:, :, 0].max()),
            ('E_MIN', ai[:, :, 1].min()),
            ('W_MAX', ai[:, :, 2].max()),
        ])
        return vals


class StatsRelevantControl(object):
    tag = "RC"

    def __init__(self):
        self.analyzer = None

    def init_lineage(self, rep, lin):
        targ = lin.targets[0]
        tset = targ.calc_distinct_outputs()
        self.analyzer = RelevantControlAnalyzer(lin.world, tset)
        self.regs = lin.params.reg_channels

    def calc_stats(self, pop):
        rc = self.analyzer.numpy_info_from_collection(pop)
        assert not np.any(np.isnan(rc.ravel()))

        # Summarize across the population
        ameans = np.mean(rc, axis=0)

        vals = []

        # Record the mean of all information measures
        regs = self.regs
        for i, c in enumerate(range(regs)):
            vals.append(('C_{}'.format(c + 1), ameans[i]))

        reg_mean = ameans.mean(axis=0)
        vals.extend([
            ('C_MEAN', reg_mean),
            ('C_MAX', rc.max()),
            ('C_MNMX', ameans.max()),
        ])
        return vals


class StatsFitness(object):
    tag = "F"

    def __init__(self):
        self.fits = None

    def init_lineage(self, rep, lin):
        self.fits = np.zeros(lin.params.population_size)

    def calc_stats(self, pop):
        pop.get_fitnesses(self.fits)
        return [
            ('MEAN', self.fits.mean()),
            ('VAR', self.fits.var()),
            ('MAX', self.fits.max()),
        ]


class StatsMutualInformation(object):
    tag = "MI"

    def __init__(self):
        self.analyzer = None
        self.regs = None

    def init_lineage(self, rep, lin):
        assert isinstance(lin, FullLineage)
        # TODO: Should really use the target that is configured in the
        # generations.
        self.target = lin.targets[0]
        self.regs = lin.params.reg_channels
        self.analyzer = MutualInfoAnalyzer(lin.world, self.target.calc_categories())

    def calc_stats(self, pop):
        mi = self.analyzer.numpy_info_from_collection(pop)

        # Reshape to remove the extra layer (there is only ONE environment)
        mi.shape = pop.size, self.regs

        # Summarize across the population
        ameans = mi.mean(axis=0)

        vals = []

        # Record the mean of all information measures
        regs = self.regs
        for i, c in enumerate(range(regs)):
            vals.append(('{}'.format(c + 1), ameans[i]))

        vals.extend([
            ('MEAN', ameans.mean()),
            ('MAX', mi.max()),
        ])
        return vals


class StatsNeighbourhood(object):
    tag = "NB"

    def __init__(self, sample_per_net=20, one_step_proportion=.5):
        self.fits = None
        self.sample_per_net = sample_per_net
        self.one_step_proportion = one_step_proportion

    def init_lineage(self, rep, lin):
        assert isinstance(lin, FullLineage)
        self.target = lin.targets[0]

    def calc_stats(self, pop):
        nayb = PopulationNeighbourhood(
            pop, self.sample_per_net, self.one_step_proportion)
        self.target.assess_collection(nayb.neighbours)
        fits = nayb.neighbours.fitnesses
        n1 = sum(fits == 1.0)
        perc = float(n1) / float(len(fits))
        return [
            ('PERC', perc),
            ('MEAN', fits.mean()),
            ('VAR', fits.var()),
            ('MED', np.median(fits)),
        ]


class StatsBindings(object):
    tag = "BD"

    def __init__(self):
        self.fits = None

    def init_lineage(self, rep, lin):
        pass

    def calc_stats(self, pop):
        bindings = pop.active_bindings
        return [
            ('MEAN', bindings.mean()),
            ('VAR', bindings.var()),
            ('MIN', bindings.min()),
        ]


class StatsRobustness(object):
    tag = "RB"

    def __init__(self):
        self.target = None

    def init_lineage(self, rep, lin):
        self.target = lin.targets[0]

    def calc_stats(self, pop):
        self.target.assess_collection(pop)
        p_fit = pop.fitnesses
        pop_mean = p_fit.mean()

        nay = PopulationNeighbourhood(pop, 100)
        nay_coll = nay.neighbours
        self.target.assess_collection(nay_coll)

        # How many
        nay_fit = nay_coll.fitnesses
        better_than_mean = np.where(nay_fit >= pop_mean)[0].size
        prop_better = float(better_than_mean) / float(nay_coll.size)
        perfect = float(np.where(pop.fitnesses == 1.0)[0].size) / float(pop.size)
        mutated_perfect = float(np.where(nay_fit == 1.0)[0].size) / float(nay_coll.size)

        return [
            ('PROP', prop_better),
            ('MEAN', pop_mean),
            ('BEST_POP', perfect),
            ('BEST_MUT', mutated_perfect),
        ]


class StatsLag(object):
    def __init__(self, experiment):
        self.experiment = experiment
        self.session = experiment.database.session
        self.replicate = None
        self.lineage = None
        self.ac_analyzer = None
        self.mi_analyzer = None
        self.done = {}

        self.first_best = 'FIRST_BEST'
        self.first_control = 'FIRST_CONTROL'

        # TODO: should only load relevant ones
        for srep in self.session.query(StatsReplicateRecord).all():
            self.done[(srep.treatment_id, srep.replicate_id, srep.kind)] = srep

    def is_done(self, rep, kind):
        return (rep.treatment.seq, rep.seq, kind) in self.done

    def visit_lineage(self, rep, lin):
        log.info("{}".format(rep)).push().add()

        # Setup
        self.replicate = rep
        self.lineage = lin
        self.ac_analyzer = AverageControlAnalyzer(lin.world)
        self.mi_analyzer = MutualInfoAnalyzer(lin.world, lin.targets[0].calc_categories())

        # Analysis
        if not self.is_done(rep, self.first_best):
            fgen = self.find_first_winner()
            self.session.add(StatsReplicateRecord(rep, self.first_best, fgen))

        if not self.is_done(rep, self.first_control):
            cgen = self.find_first_control()
            self.session.add(StatsReplicateRecord(rep, self.first_control, cgen))

        self.session.commit()

        log.pop()

    def find_first_winner(self):
        first_win = None
        for g in self.lineage._generations.where("best == 1.0"):
            first_win = g['generation']
            break

        log.info("Found first winner at generation {}".format(first_win))
        return first_win

    def find_first_control(self):
        # Load last generation and check if we got a master gene
        c_index = self.find_controlled(self.lineage.population)
        if c_index is None:
            log.debug("No control evolved")
            return None
        log.debug("We got a master gene at the last generation.")

        # Grab the first one and load the ancestry
        net = self.lineage.population[c_index]
        first_control = self.get_ancestry_control(net)

        # Just one
        log.info("First control is at generation {}.".format(first_control))
        return first_control

    def find_controlled(self, collection):
        output_size = self.lineage.params.out_channels

        mutual = self.mi_analyzer.numpy_info_from_collection(collection)
        # Ditch the empty dimension
        mutual.shape = mutual.shape[:-1]

        ai = self.ac_analyzer.numpy_info_from_collection(collection)
        control = ai[:, :, :output_size]
        entropy = ai[:, :, output_size:]

        for i, (m_net, c_net, e_net) in enumerate(zip(mutual, control, entropy)):
            for j, (m_reg, c_reg, e_reg) in enumerate(zip(m_net, c_net, e_net)):
                # It must have full information about the environment
                if m_reg == 1.0:
                    # And this must translate into control
                    if (c_reg == 1.0).all():
                            # The control must explain all of the entropy
                            if (e_reg == c_reg).all():
                                return i
        return None

    def get_ancestry_control(self, network):
        log.debug("loading ancestry of winner {}".format(network.identifier))
        anc = self.lineage.get_ancestry(network.identifier)

        log.debug("Calculating ancestry control of {} networks.".format(anc.size))
        c_index = self.find_controlled(anc)

        # We know that at least the last one should be good!
        assert c_index is not None

        first_ancestor = anc[c_index]

        # Make sure the network has a fitnes
        self.lineage.targets[0].assess(first_ancestor)
        log.debug("First ancestor with master gene is at {}.".format(
            first_ancestor.generation))
        self.replicate.draw_net('first-control', first_ancestor, target=self.lineage.targets[0])
        return first_ancestor.generation


class StatsGenerations(object):
    def __init__(self, experiment):
        self.experiment = experiment
        self.session = experiment.database.session
        self.max_gen = 'MAX_GENERATIONS'
        self.done = {}

        # TODO: should only load relevant ones
        for srep in self.session.query(StatsReplicateRecord).all():
            self.done[(srep.treatment_id, srep.replicate_id, srep.kind)] = srep

    def is_done(self, rep, kind):
        return (rep.treatment.seq, rep.seq, kind) in self.done

    def visit_lineage(self, rep, lin):
        log.info("{}".format(rep)).push().add()

        # Analysis
        if not self.is_done(rep, self.max_gen):
            fgen = lin.generation
            self.session.add(StatsReplicateRecord(rep, self.max_gen, fgen))
            log.info("Final generation is {}".format(fgen))

        self.session.commit()
        log.pop()
