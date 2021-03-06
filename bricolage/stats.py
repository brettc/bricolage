"""High level analysis tools
"""
import logtools
import numpy as np
from .analysis_ext import (
    MutualInfoAnalyzer,
    OutputControlAnalyzer,
    RelevantControlAnalyzer,
    FastCAnalyzer,
    MIAnalyzer,
    WCAnalyzer,
    FastCandBAnalyzer,
)
from .analysis import AverageControlAnalyzer
from .lineage import FullLineage
from .experimentdb import (
    StatsGroupRecord,
    StatsRecord,
    StatsReplicateRecord,
    NetworkRecord,
)
from .neighbourhood import PopulationNeighbourhood
from bricolage.experiment import Experiment
from bricolage.graph_maker import SignalFlowGraph
from bricolage.analysis_ext import NetworkAnalysis
import inspect
from logic2 import modules_changed

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
            if inspect.isclass(kls):
                k = kls()
            else:
                k = kls
            # We need a tag attribute.
            assert hasattr(k, "tag")
            self.todo.add(k)

        # Load all stats groups already done.
        # NOTE: Not sure this is the best way to do it yet...
        self.done = {}

        log.info("Loading all stats groups...")
        for grp in self.session.query(StatsGroupRecord).all():
            self.done[
                (grp.treatment_id, grp.replicate_id, grp.generation, grp.tag)
            ] = grp
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
        log.info(
            "Doing stats for generation {}, {}".format(self.replicate.path, gen_num)
        )
        for stats in self.todo:
            if self.is_done(self.replicate, gen_num, stats.tag):
                log.info(
                    "Skipping already created group for %s, %s", gen_num, stats.tag
                )
                continue

            self.add_stats_group(self.replicate, gen_num, stats.tag)
            named_values = stats.calc_stats(pop)
            self.session.add_all(
                [
                    StatsRecord(
                        self.replicate,
                        gen_num,
                        "{}_{}".format(stats.tag, n),
                        v,
                        stats.tag,
                    )
                    for n, v in named_values
                ]
            )

        self.session.commit()


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
                vals.append(("C_{}_{}".format(c + 1, j + 1), control[i, j]))
                vals.append(("E_{}_{}".format(c + 1, j + 1), entropy[i, j]))

        vals.extend([("C_MEAN", control.mean()), ("E_MEAN", entropy.mean())])
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
            vals.append(("C_{}".format(c + 1), ameans[i, 0]))
            vals.append(("E_{}".format(c + 1), ameans[i, 1]))
            vals.append(("W_{}".format(c + 1), ameans[i, 2]))

        reg_mean = ameans.mean(axis=0)
        vals.extend(
            [
                ("C_MEAN", reg_mean[0]),
                ("E_MEAN", reg_mean[1]),
                ("W_MEAN", reg_mean[2]),
                ("C_MAX", ai[:, :, 0].max()),
                ("E_MIN", ai[:, :, 1].min()),
                ("W_MAX", ai[:, :, 2].max()),
            ]
        )
        return vals


class StatsRelevantControl(object):
    def __init__(self, use_natural=True):
        if use_natural:
            self.tag = "AC"
        else:
            self.tag = "PC"

        self.analyzer = None
        self.use_natural = use_natural

    def init_lineage(self, rep, lin):
        targ = lin.targets[0]
        tset = targ.calc_distinct_outputs()
        self.analyzer = RelevantControlAnalyzer(
            lin.world, tset, use_natural=self.use_natural
        )
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
            vals.append(("{}".format(c + 1), ameans[i]))

        reg_mean = ameans.mean(axis=0)
        vals.extend([("MEAN", reg_mean), ("MAX", rc.max()), ("MNMX", ameans.max())])
        return vals


class StatsBaseControl(object):
    def __init__(self, tag, indexes, target1, target2):
        self.tag = tag
        self.indexes = indexes
        self.target1 = target1
        self.target2 = target2
        self.analyzer = None

    def init_lineage(self, rep, lin):
        self.analyzer = FastCAnalyzer(
            lin.world, self.indexes, self.target1, self.target2
        )
        self.regs = lin.params.reg_channels

    def calc_stats(self, pop):
        rc = self.analyzer.analyse_collection(pop)
        assert not np.any(np.isnan(rc.ravel()))

        # Summarize across the population
        ameans = np.mean(rc, axis=0)

        vals = []

        # Record the mean of all information measures
        regs = self.regs
        for i, c in enumerate(range(regs)):
            vals.append(("{}".format(c + 1), ameans[i]))

        reg_mean = ameans.mean(axis=0)
        vals.extend([("MEAN", reg_mean), ("MAX", rc.max()), ("MNMX", ameans.max())])
        return vals


class StatsFastControl(StatsBaseControl):
    def __init__(self, tag, indexes, target1, target2):
        super(StatsFastControl, self).__init__(tag, indexes, target1, target2)

    def init_lineage(self, rep, lin):
        self.analyzer = FastCAnalyzer(
            lin.world, self.indexes, self.target1, self.target2
        )
        self.regs = lin.params.reg_channels


class StatsFastWithBackgroundControl(StatsBaseControl):
    def __init__(self, tag, indexes, target1, target2):
        super(StatsFastWithBackgroundControl, self).__init__(
            tag, indexes, target1, target2
        )

    def init_lineage(self, rep, lin):
        self.analyzer = FastCandBAnalyzer(
            lin.world, self.indexes, self.target1, self.target2
        )
        self.regs = lin.params.reg_channels


class StatsWeightedControl(object):
    def __init__(self, tag, indexes, target1, target2, weighting):
        self.tag = tag
        self.indexes = indexes
        self.target1 = target1
        self.target2 = target2
        self.weighting = weighting
        self.analyzer = None

    def init_lineage(self, rep, lin):
        self.analyzer = WCAnalyzer(
            lin.world, self.indexes, self.target1, self.target2, self.weighting
        )
        self.regs = lin.params.reg_channels

    def calc_stats(self, pop):
        rc = self.analyzer.analyse_collection(pop)
        assert not np.any(np.isnan(rc.ravel()))

        # Summarize across the population
        ameans = np.mean(rc, axis=0)

        vals = []

        # Record the mean of all information measures
        regs = self.regs
        for i, c in enumerate(range(regs)):
            vals.append(("{}".format(c + 1), ameans[i]))

        reg_mean = ameans.mean(axis=0)
        vals.extend([("MEAN", reg_mean), ("MAX", rc.max()), ("MNMX", ameans.max())])
        return vals


class StatsFitness(object):
    tag = "F"

    def __init__(self):
        self.fits = None

    def init_lineage(self, rep, lin):
        self.fits = np.zeros(lin.params.population_size)

    def calc_stats(self, pop):
        fits = np.asarray(pop.fitnesses)
        return [("MEAN", fits.mean()), ("VAR", fits.var()), ("MAX", fits.max())]


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
            vals.append(("{}".format(c + 1), ameans[i]))

        vals.extend(
            [("MEAN", ameans.mean()), ("MAX", mi.max()), ("MNMX", ameans.max())]
        )
        return vals


class StatsNeighbourhood(object):
    tag = "NB"

    def __init__(self, sample_per_net=20, one_step_proportion=0.5):
        self.fits = None
        self.sample_per_net = sample_per_net
        self.one_step_proportion = one_step_proportion

    def init_lineage(self, rep, lin):
        assert isinstance(lin, FullLineage)
        self.target = lin.targets[0]

    def calc_stats(self, pop):
        nayb = PopulationNeighbourhood(
            pop, self.sample_per_net, self.one_step_proportion
        )
        self.target.assess_collection(nayb.neighbours)
        fits = nayb.neighbours.fitnesses
        n1 = sum(fits == 1.0)
        perc = float(n1) / float(len(fits))
        return [
            ("PERC", perc),
            ("MEAN", fits.mean()),
            ("VAR", fits.var()),
            ("MED", np.median(fits)),
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
            ("MEAN", bindings.mean()),
            ("VAR", bindings.var()),
            ("MIN", bindings.min()),
        ]


class StatsRobustness(object):
    tag = "RB"

    def __init__(self):
        self.target = None

    def init_lineage(self, rep, lin):
        self.target = lin.targets[0]

    def calc_stats(self, pop):
        self.target.assess_collection(pop)
        p_fit = np.asarray(pop.fitnesses)
        pop_mean = p_fit.mean()

        nay = PopulationNeighbourhood(pop, 100)
        nay_coll = nay.neighbours
        nay_fit = np.asarray(self.target.assess_collection(nay_coll))
        better_than_mean = np.where(nay_fit >= pop_mean)[0].size
        prop_better = float(better_than_mean) / float(nay_coll.size)
        perfect = float(np.where(pop.fitnesses == 1.0)[0].size) / float(pop.size)
        mutated_perfect = float(np.where(nay_fit == 1.0)[0].size) / float(nay_coll.size)

        return [
            ("PROP", prop_better),
            ("MEAN", pop_mean),
            ("BEST_POP", perfect),
            ("BEST_MUT", mutated_perfect),
        ]


class StatsEnvironmental(object):
    tag = "EB"

    def __init__(self):
        self.target = None

    def init_lineage(self, rep, lin):
        self.target = lin.targets[0]

    def calc_stats(self, pop):
        self.target.assess_collection(pop)
        p_fit = np.asarray(pop.fitnesses)

        robusts = pop.robustness()
        rob_mean = robusts.mean()

        rob_scaled = p_fit * robusts
        robs_mean = rob_scaled.mean()

        return [("MEAN", rob_mean), ("FIT_MEAN", robs_mean)]


class StatsLag(object):
    def __init__(self, experiment):
        self.experiment = experiment
        self.session = experiment.database.session
        self.replicate = None
        self.lineage = None
        self.rc_analyzer = None
        self.mi_analyzer = None
        self.done = {}

        self.first_best = "FIRST_BEST"
        self.first_control = "FIRST_CONTROL"
        self.first_master = "FIRST_MASTER"

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
        targ = lin.targets[0]
        tset = targ.calc_distinct_outputs()
        self.rc_analyzer = RelevantControlAnalyzer(lin.world, tset)
        self.mi_analyzer = MutualInfoAnalyzer(lin.world, targ.calc_categories())

        # Analysis
        if not self.is_done(rep, self.first_best):
            fgen = self.find_first_winner()
            self.session.add(StatsReplicateRecord(rep, self.first_best, fgen))

        if not self.is_done(rep, self.first_control):
            cgen, mgen = self.find_first_control()
            self.session.add(StatsReplicateRecord(rep, self.first_control, cgen))
            self.session.add(StatsReplicateRecord(rep, self.first_master, mgen))

        self.session.commit()

        log.pop()

    def find_first_winner(self):
        first_win = None
        for g in self.lineage._generations.where("best == 1.0"):
            first_win = g["generation"]
            break

        log.info("Found first winner at generation {}".format(first_win))
        return first_win

    def find_first_control(self):
        # Load last generation and check if we got a master gene
        _, m_index = self.find_controlled(self.lineage.population)
        if m_index is None:
            log.info("No master evolved")
            return None, None

        log.info("We got a master gene at the last generation.")

        # Grab the first one and load the ancestry
        net = self.lineage.population[m_index]
        first_control, first_master = self.get_ancestry_control(net)

        # Just one
        log.info("First control is at generation {}.".format(first_master))
        return first_control, first_master

    def find_controlled(self, collection):
        """Control requires both be 1.0 too"""
        rc = self.rc_analyzer.numpy_info_from_collection(collection)
        rc_indexes = np.where(rc == 1.0)[0]
        if rc_indexes.size > 0:
            first_control = rc_indexes[0]
        else:
            first_control = None

        # Yes, there is a faster way...
        mi = self.mi_analyzer.numpy_info_from_collection(collection)

        # TODO: fix this stupidity
        mi.shape = mi.shape[:2]
        mi_and_rc = mi * rc
        master_indexes = np.where(mi_and_rc == 1.0)[0]
        if master_indexes.size > 0:
            first_master = master_indexes[0]
        else:
            first_master = None

        return first_control, first_master

    def get_ancestry_control(self, network):
        log.debug("loading ancestry of winner {}".format(network.identifier))
        anc = self.lineage.get_ancestry(network.identifier)

        log.debug("Calculating ancestry control of {} networks.".format(anc.size))
        control, master = self.find_controlled(anc)

        # We know that at least the last one should be good!
        assert master is not None

        first_control = anc[control]
        first_master = anc[master]

        # Make sure the network has a fitness
        log.debug(
            "First ancestor with master gene is at {}.".format(first_master.generation)
        )
        self.replicate.draw_net(
            "first-control", first_master, target=self.lineage.targets[0]
        )
        return first_control.generation, first_master.generation


class StatsGenerations(object):
    def __init__(self, experiment):
        self.experiment = experiment
        self.session = experiment.database.session
        self.max_gen = "MAX_GENERATIONS"
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


class StatsTime(object):
    def __init__(self, experiment, streak_length):
        self.experiment = experiment
        self.session = experiment.database.session
        self.replicate = None
        self.lineage = None
        self.done = {}
        self.streak_length = streak_length

        self.last_generation = "LAST_GENERATION"
        self.begin_streak = "BEGIN_STREAK"

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

        # Analysis
        if not self.is_done(rep, self.last_generation):
            fgen = lin.generation
            log.info("Last generation is {}".format(fgen))
            self.session.add(StatsReplicateRecord(rep, self.last_generation, fgen))

        if not self.is_done(rep, self.begin_streak):
            l = lin.final_streak_length(self.streak_length + 1)
            beg = lin.generation - l + 1
            log.info("Last streak begins at generation {}".format(beg))
            self.session.add(StatsReplicateRecord(rep, self.begin_streak, beg))

        self.session.commit()
        log.pop()


class StatsMI(object):
    def __init__(self, tag, indexes, target_num=0):
        self.tag = tag
        self.target_num = target_num
        self.indexes = indexes
        self.analyzer = None
        self.regs = None

    def init_lineage(self, rep, lin):
        assert isinstance(lin, FullLineage)

        self.regs = lin.params.reg_channels
        target = lin.targets[self.target_num]
        categories = target.calc_categories(self.indexes)
        assert set(categories) == set([0, 1])
        self.analyzer = MIAnalyzer(lin.world, categories)

    def calc_stats(self, pop):
        mi = self.analyzer.analyse_collection(pop)

        # Summarize across the population
        ameans = mi.mean(axis=0)

        vals = []

        # Record the mean of all information measures
        regs = self.regs
        for i, c in enumerate(range(regs)):
            vals.append(("{}".format(c + 1), ameans[i]))

        vals.extend(
            [("MEAN", ameans.mean()), ("MAX", mi.max()), ("MNMX", ameans.max())]
        )
        return vals


class StatsMaster(object):
    """Work out how many master genes there are"""

    def __init__(self, tag, indexes, target1, target2, target_num=0, bowties_too=True):
        self.tag = tag
        self.indexes = indexes
        self.target1 = target1
        self.target2 = target2
        self.target_num = target_num
        self.target = None
        self.r_analyzer = None
        self.m_analyzer = None
        self.bowties_too = bowties_too

    def init_lineage(self, rep, lin):
        self.r_analyzer = FastCAnalyzer(
            lin.world, self.indexes, self.target1, self.target2
        )
        self.target = lin.targets[self.target_num]
        categories = self.target.calc_categories(self.indexes)
        self.m_analyzer = MIAnalyzer(lin.world, categories)
        self.regs = lin.params.reg_channels

    def do_bowties(self, pop):
        log.info("Calculating Bowties...")
        bowties = np.zeros(pop.size, np.bool_)
        fit_bowties = np.zeros(pop.size, np.bool_)

        for i, net in enumerate(pop):
            ana = NetworkAnalysis(net)
            fg = SignalFlowGraph(ana)
            if fg.has_bowtie():
                bowties[i] = 1
                if net.fitness == 1.0:
                    fit_bowties[i] = 1

        return bowties, fit_bowties

    def calc_stats(self, pop):
        # Get the fitnesses
        self.target.assess_collection(pop)
        fits = pop.fitnesses

        # Analyse
        log.info("Calculating Information...")
        rc = self.r_analyzer.analyse_collection(pop)
        mi = self.m_analyzer.analyse_collection(pop)

        def any_true_for_net(bool_arr):
            # Is anything true?
            return bool_arr.sum(axis=1) >= 1

        def freq(bool_arr):
            return bool_arr.sum() / float(pop.size)

        def gene_freq(bool_arr):
            return bool_arr.sum() / float(rc.size)

        gene_info = mi == 1.0
        gene_control = rc == 1.0
        gene_master = gene_info & gene_control

        # Look at what is going on in combinations
        nets_fit = fits == 1.0
        nets_with_info = any_true_for_net(gene_info)
        nets_with_control = any_true_for_net(gene_control)
        nets_with_master = any_true_for_net(gene_info & gene_control)
        nets_fit_master = nets_fit & nets_with_master
        nets_fit_info = nets_with_info & nets_fit
        nets_fit_control = nets_with_control & nets_fit

        vals = [
            ("FIT", freq(nets_fit)),
            ("INFO", freq(nets_with_info)),
            ("CONTROL", freq(nets_with_control)),
            ("MASTER", freq(nets_with_master)),
            ("FIT_MASTER", freq(nets_fit_master)),
            ("FIT_INFO", freq(nets_fit_info)),
            ("FIT_CONTROL", freq(nets_fit_control)),
            ("GENE_INFO", gene_freq(gene_info)),
            ("GENE_CONTROL", gene_freq(gene_control)),
            ("GENE_MASTER", gene_freq(gene_master)),
        ]

        if self.bowties_too:

            bowties, fit_bowties = self.do_bowties(pop)

            master_no_bow = nets_with_master & ~bowties
            fit_master_no_bow = nets_fit_master & ~bowties
            bowtie_no_master = bowties & ~nets_with_master
            fit_bowtie_no_master = fit_bowties & ~nets_with_master
            vals.extend(
                [
                    ("BOW_PROB", freq(bowties)),
                    ("BOW_FIT_PROB", freq(fit_bowties)),
                    ("MASTER_NO_BOW", freq(master_no_bow)),
                    ("FIT_MASTER_NO_BOW", freq(fit_master_no_bow)),
                    ("BOW_NO_MASTER", freq(bowtie_no_master)),
                    ("FIT_BOW_NO_MASTER", freq(fit_bowtie_no_master)),
                ]
            )

        return vals


class StatsRCMaster(object):
    """Work out how many master genes there are"""

    tag = "RCM"

    def __init__(self, target_num=0):
        self.target_num = target_num
        self.target = None
        self.r_analyzer = None
        self.m_analyzer = None

    def init_lineage(self, rep, lin):
        self.target = lin.targets[self.target_num]
        tset = self.target.calc_distinct_outputs()
        self.r_analyzer = RelevantControlAnalyzer(lin.world, tset)
        categories = self.target.calc_categories()
        self.m_analyzer = MutualInfoAnalyzer(lin.world, categories)
        self.regs = lin.params.reg_channels

    def calc_stats(self, pop):
        # Get the fitnesses
        self.target.assess_collection(pop)
        fits = pop.fitnesses

        # Analyse
        rc = self.r_analyzer.numpy_info_from_collection(pop)
        mi = self.m_analyzer.numpy_info_from_collection(pop)
        mi.shape = mi.shape[:-1]

        def any_true_for_net(bool_arr):
            # Is anything true?
            return bool_arr.sum(axis=1) >= 1

        def freq(bool_arr):
            return bool_arr.sum() / float(pop.size)

        def gene_freq(bool_arr):
            return bool_arr.sum() / float(rc.size)

        gene_info = mi == 1.0
        gene_control = rc == 1.0
        gene_master = gene_info & gene_control

        # Look at what is going on in combinations
        nets_fit = fits == 1.0
        nets_with_info = any_true_for_net(gene_info)
        nets_with_control = any_true_for_net(gene_control)
        nets_with_master = any_true_for_net(gene_info & gene_control)
        nets_fit_master = nets_fit & nets_with_master
        nets_fit_info = nets_with_info & nets_fit
        nets_fit_control = nets_with_control & nets_fit

        vals = [
            ("FIT", freq(nets_fit)),
            ("INFO", freq(nets_with_info)),
            ("CONTROL", freq(nets_with_control)),
            ("MASTER", freq(nets_with_master)),
            ("FIT_MASTER", freq(nets_fit_master)),
            ("FIT_INFO", freq(nets_fit_info)),
            ("FIT_CONTROL", freq(nets_fit_control)),
            ("GENE_INFO", gene_freq(gene_info)),
            ("GENE_CONTROL", gene_freq(gene_control)),
            ("GENE_MASTER", gene_freq(gene_master)),
        ]
        return vals


class StatsCisInDegree(object):
    tag = "CIS"

    def __init__(self):
        pass

    def init_lineage(self, rep, lin):
        pass

    def calc_stats(self, pop):
        cis = pop.active_cis()
        return [("{:02d}".format(i), n) for (i, n) in enumerate(cis)]


class StatsBowtie(object):
    tag = "BOW"

    def __init__(self, target_num=0):
        self.target_num = target_num

    def init_lineage(self, rep, lin):
        self.target = lin.targets[self.target_num]

    def calc_stats(self, pop):
        # Get the fitnesses
        self.target.assess_collection(pop)

        bowties = 0
        fit_bowties = 0
        for net in pop:
            ana = NetworkAnalysis(net)
            fg = SignalFlowGraph(ana)
            if fg.has_bowtie():
                bowties += 1
                if net.fitness == 1.0:
                    fit_bowties += 1

        log.info("Successful Bowties: {} - fit {}".format(bowties, fit_bowties))
        return [
            ("PROB", bowties / float(pop.size)),
            ("FIT_PROB", fit_bowties / float(pop.size)),
        ]


class StatsFirstMaster(object):
    """Work out how many master genes there are"""

    def __init__(
        self, experiment, tag, indexes, target1, target2, weighting, target_num=0
    ):
        self.session = experiment.database.session
        self.tag = tag
        self.indexes = indexes
        self.target1 = target1
        self.target2 = target2
        self.weighting = weighting
        self.target_num = target_num

        self.first_master = self.tag + "_FIRST_MASTER"
        self.r_analyzer = None
        self.m_analyzer = None
        self.replicate = None
        self.done = {}

        # TODO: should only load relevant ones
        for srep in self.session.query(StatsReplicateRecord).all():
            self.done[(srep.treatment_id, srep.replicate_id, srep.kind)] = srep

    def is_done(self, rep, kind):
        return (rep.treatment.seq, rep.seq, kind) in self.done

    def visit_lineage(self, rep, lin):
        self.replicate = rep
        log.info("{}".format(rep)).push().add()
        if not self.is_done(rep, self.first_master):
            fgen = self.find_first_master(lin)
            self.session.add(StatsReplicateRecord(rep, self.first_master, fgen))
            self.session.commit()
        log.pop()

    def find_first_master(self, lin):
        self.r_analyzer = WCAnalyzer(
            lin.world, self.indexes, self.target1, self.target2, self.weighting
        )
        target = lin.targets[self.target_num]
        categories = target.calc_categories(self.indexes)
        self.m_analyzer = MIAnalyzer(lin.world, categories)

        # Is there any in the last population?
        last_pop_i = self.master_indexes(lin.population)
        if len(last_pop_i) == 0:
            log.info("No Master")
            return None

        log.info("Found master in final pop, loading ancestry")

        net = lin.population[last_pop_i[0]]

        # Load the lineage
        anc = lin.get_ancestry(net.identifier)
        log.info("Ancestry size is {}".format(anc.size))
        anc_i = self.master_indexes(anc)
        first_master = anc[anc_i[0]]
        log.info("Master at {}".format(first_master.generation))
        target.assess_collection(anc)
        self.replicate.draw_net("first-control", first_master, target=target)
        return first_master.generation

    def master_indexes(self, coll):
        rc = self.r_analyzer.analyse_collection(coll)
        mi = self.m_analyzer.analyse_collection(coll)

        res = ((rc == 1.0) & (mi == 1.0)).sum(axis=1)
        return np.where(res >= 1)[0]


class StatsNetworkDiffs(object):
    tag = "MUT"

    def __init__(self, experiment):
        self.experiment = experiment
        self.session = experiment.database.session
        self.replicate = None
        self.lineage = None
        self.done = {}

        self.gens = "GEN_DIFF"
        self.edges = "EDGE_DIFF"

        # TODO: should only load relevant ones
        for srep in self.session.query(StatsReplicateRecord).all():
            self.done[(srep.treatment_id, srep.replicate_id, srep.kind)] = srep

    def is_done(self, rep, kind):
        return (rep.treatment.seq, rep.seq, kind) in self.done

    def visit_lineage(self, rep, lin):
        log.info("{}".format(rep))

        # # Analysis
        if self.is_done(rep, self.gens):
            return

        gfirst = lin.first_winning_generation()
        log.info("First generation is {}".format(gfirst))
        if gfirst is not None:

            gen = lin.get_generation(gfirst)
            log.info("First generation is {}".format(gen))
            cur_net = gen.get_best(1)[0]

            gen_diff = cur_net.generation
            log.info("Generations for change {}".format(gen_diff))

            anc = lin.get_ancestry(cur_net.identifier)
            first_net = anc[0]
            log.info("Loaded ancestry")

            cur_ana = NetworkAnalysis(cur_net)
            first_ana = NetworkAnalysis(first_net)

            cur_edges = cur_ana.get_active_edges()
            first_edges = first_ana.get_active_edges()

            edge_diff = len(cur_edges ^ first_edges)
            log.info("Causal diffs {}".format(edge_diff))
        else:
            gen_diff = None
            edge_diff = None

        self.session.add(StatsReplicateRecord(rep, self.gens, gen_diff))
        self.session.add(StatsReplicateRecord(rep, self.edges, edge_diff))
        self.session.commit()


class StatsNetworkChanges(object):
    kind = "MC"

    def __init__(self, experiment):
        self.experiment = experiment
        self.session = experiment.database.session
        self.done = set()

        for srep in (
            self.session.query(NetworkRecord)
            .filter(NetworkRecord.kind == self.kind)
            .all()
        ):
            self.done.add((srep.treatment_id, srep.replicate_id))

    def wants_replicate(self, rep):
        wanted = (rep.treatment.seq, rep.seq) not in self.done
        return wanted

    def visit_lineage(self, rep, lin):
        log.info("{}".format(rep))

        glast = lin.generation
        log.info("Last generation is {}".format(glast))

        cur_net = lin.population.get_best(2)[0]
        if cur_net.fitness == 1.0:
            anc = lin.get_ancestry(cur_net.identifier)
            first_net = anc[0]
            changes = modules_changed(first_net, cur_net)
            log.info("Total changes = {}".format(len(changes)))

            for g, m in changes:
                rec = NetworkRecord(rep, cur_net, g, m, self.kind, 1.0)
                self.session.add(rec)
        else:
            log.info("Didn't find a successful network")
            # Just add a "failure" record
            rec = NetworkRecord(rep, cur_net, -1, -1, self.kind, 0.0)

        self.session.commit()


class StatsLastWinner(object):
    def __init__(self, experiment):
        self.experiment = experiment
        self.session = experiment.database.session
        self.done = set()
        self.kind = "LAST_WINNER"

        for srep in (
            self.session.query(StatsReplicateRecord)
            .filter(StatsReplicateRecord.kind == self.kind)
            .all()
        ):
            self.done.add((srep.treatment_id, srep.replicate_id))

    def wants_replicate(self, rep):
        wanted = (rep.treatment.seq, rep.seq) not in self.done
        return wanted

    def visit_lineage(self, rep, lin):
        log.info("{}".format(rep)).push().add()

        w, b = lin.population.worst_and_best()
        if b == 1.0:
            fgen = lin.generation
        else:
            fgen = None

        log.info("Last generation is {}".format(fgen))
        self.session.add(StatsReplicateRecord(rep, self.kind, fgen))
        self.session.commit()
        log.pop()


class StatsSwitchboard(object):
    """Work out how many master genes there are"""

    def __init__(self, tag, sample=0, target_num=0):
        self.tag = tag
        self.target_num = target_num
        self.target = None
        self.r_analyzer = None
        self.m_analyzer = None
        self.sample = sample

    def init_lineage(self, rep, lin):
        self.target = lin.targets[self.target_num]
        tset = self.target.calc_distinct_outputs()
        self.r_analyzer = RelevantControlAnalyzer(lin.world, tset)
        categories = self.target.calc_categories()
        self.m_analyzer = MutualInfoAnalyzer(lin.world, categories)
        self.regs = lin.params.reg_channels

    def calc_stats(self, pop):
        # Do we want to extend it?
        if self.sample <= 0:
            coll = pop
        else:
            pn = PopulationNeighbourhood(pop, self.sample)
            coll = pn.neighbours

        # Calculate the information
        rc = self.r_analyzer.numpy_info_from_collection(coll)
        mi = self.m_analyzer.numpy_info_from_collection(coll)

        # Get upstream info, and downstream information
        # TODO: get rid of ugliness.
        mi.shape = mi.shape[:-1]

        def mean_freq_is_1(arr):
            freqs = (arr == 1.0).sum(axis=1) / float(self.regs)
            return freqs.mean()

        upstream = mean_freq_is_1(mi)
        downstream = mean_freq_is_1(rc)
        both = upstream * downstream
        log.info("Up: {}, Down: {}, Both: {}".format(upstream, downstream, both))

        vals = [("UPSTREAM", upstream), ("DOWNSTREAM", downstream), ("BOTH", both)]
        return vals
