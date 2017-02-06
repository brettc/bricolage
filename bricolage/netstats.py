"""Stats visitors for looking around networks
"""
import logtools
# import pandas as pd
import numpy as np
from .logic2 import modules_changed
from .analysis_ext import RelevantControlAnalyzer, MutualInfoAnalyzer, AverageControlAnalyzer
from .experimentdb import GeneMeasureRecord
from .neighbourhood import NetworkNeighbourhood
import random
from enum import IntEnum


log = logtools.get_logger()


class MutationType(IntEnum):
    SAME = 0
    NOVEL = 1
    DEAD = 2

# def is_novel(net, pats, targ):
#     env_len, pat_len = targ.shape
#     if np.allclose(net.rates, targ):
#         return 0
#     matches = np.zeros(env_len, np.bool)
#     for p in pats:
#         matches |= (net.rates == p).sum(axis=1) == pat_len
#     if np.alltrue(matches):
#         return 1
#     return -1

def is_novel(net, pats, targ):
    env_len, pat_len = targ.shape
    if np.allclose(net.rates, targ):
        return 0
    cur_matches = np.zeros(env_len, np.bool)
    all_matches = np.zeros(env_len, np.bool)
    for p in pats:
        cur_matches = (net.rates == p).sum(axis=1) == pat_len
        # if it matches ALL of them, then there no plasticity
        if cur_matches.sum() == env_len:
            return -1
        all_matches |= cur_matches
    if np.alltrue(all_matches):
        return 1
    return -1


def get_pattern_ident(net, pats):
    ident = np.zeros(net.rates.shape[0], int)
    patlen = pats.shape[1]
    for i, p in enumerate(pats):
        found = np.where((net.rates == p).sum(axis=1) == patlen)[0]
        ident[found] = i + 1
    assert not np.any(ident == 0)
    return tuple(ident)


def get_novel(center, sample, pats, targ, one_step=1.0):
    nay = NetworkNeighbourhood(center=center, 
                               sample_size=sample, 
                               one_step_proportion=one_step)
    same = []
    nov = []
    dead = []
    for net in nay.neighbours:
        kind = is_novel(net, pats, targ)
        if kind == 0:
            same.append(net)
        elif kind == 1:
            nov.append(net)
        else:
            dead.append(net)
    return same, nov, dead


def count_novel(nov, pats):
    uniq = set()
    for net in nov:
        uniq.add(get_pattern_ident(net, pats))
    return len(uniq)


def count_all_novel(best, sample, pats, targ):
    uniq = set()
    nsame = 0
    nnov = 0
    ndead = 0
    for net in best:
        same, nov, dead = get_novel(net, sample, pats, targ)
        nsame += len(same)
        nnov += len(nov)
        ndead += len(dead)
        for nov_net in nov:
            uniq.add(get_pattern_ident(nov_net, pats))
    return nsame, nnov, ndead, len(uniq)


def count_genes(nn, nov):
    d = {}
    tot = 0
    for cur_net in nov:    
        ch = modules_changed(nn, cur_net)
        for gene, mod in ch:
            d.setdefault(gene, 0)
            d[gene] += 1
            tot += 1
    return d, tot


def append_count_genes(nn, nov, d):
    for cur_net in nov:    
        ch = modules_changed(nn, cur_net)
        # if len(ch) > 1:
        #     raise RuntimeError
        for gene, mod in ch:
            n = d.setdefault(gene, 0)
            d[gene] = n + 1
    return d


class GenePotential(object):
    def __init__(self, experiment, tag, net_samples, explore_samples, use_natural=False):
        self.experiment = experiment
        self.session = experiment.database.session
        self.tag = tag
        self.net_samples = net_samples
        self.explore_samples = explore_samples
        self.use_natural = use_natural

    # def wants_replicate(self, rep):
    #     wanted = (rep.treatment.seq, rep.seq) not in self.done
    #     return wanted
    #
    
    def init(self, lin):
        self.lineage = lin
        self.target = lin.targets[0]
        tset = self.target.calc_distinct_outputs()
        self.valid_patterns = np.asarray(list(tset))
        self.optimal_rates = self.target.as_array()
        log.info("Calculating information measure")

        # self.rz = RelevantControlAnalyzer(lin.world, tset, use_natural=self.use_natural)
        # self.rz_full = self.rz.numpy_info_from_collection(lin.population)
        # self.rz_info = self.rz_full.mean(axis=0) 

        self.rz = AverageControlAnalyzer(lin.world)
        info = np.asarray(self.rz.analyse_collection(lin.population))
        info_dim = info.shape[2] / 3
        self.rz_full = info[:,:,:info_dim].sum(axis=2)
        self.rz_info = self.rz_full.mean(axis=0)

        categories = self.target.calc_categories()
        self.mz = MutualInfoAnalyzer(lin.world, categories)
        self.mz_full = self.mz.numpy_info_from_collection(lin.population)
        self.mz_full.shape = self.mz_full.shape[:-1]
        self.mz_info = self.mz_full.mean(axis=0) 

        self.regs = lin.params.reg_channels

    def get_networks(self):
        # Just work on the last population
        best = self.lineage.population.get_best()
        uniq_idents = set()
        uniq_best = []
        for net in best:
            if net.identifier not in uniq_idents:
                uniq_idents.add(net.identifier)
                uniq_best.append(net)

        assert len(uniq_best) >= self.net_samples
        return random.sample(uniq_best, self.net_samples)

    def visit_lineage(self, rep, lin):
        log.info("{}".format(rep)).push().add()
        self.init(lin)

        nets = self.get_networks()
        log.info("Generating {} mutants for {} networks".format(self.explore_samples,
                                                                self.net_samples))
        # self.analyse_multiple(rep, nets)
        self.analyse_single(rep, nets)
        self.session.commit()
        log.pop()

    def analyse_single(self, rep, nets):
        counts = [{}, {}, {}]
        # uniq = set()

        # Assemble everything
        for cur_net in nets:
            mutants = get_novel(
                cur_net, self.explore_samples, self.valid_patterns, self.optimal_rates)
            for mut, cnt in zip(mutants, counts): 
                append_count_genes(cur_net, mut, cnt)
            # for nov_net in mutants[MutationType.NOVEL]:
            #     uniq.add(get_pattern_ident(nov_net, self.valid_patterns))

        # Produce a single data-point for each gene / replicate
        for i in range(self.regs):
            rz = self.rz_info[i]
            mz = self.mz_info[i]
            for mtype in MutationType:
                gene_counts = counts[mtype]
                count = gene_counts.get(i, 0)
                gm = GeneMeasureRecord(rep, cur_net, i + 1, self.tag, mtype.name)
                gm.measure = rz
                gm.found = count
                self.session.add(gm)

                gmm = GeneMeasureRecord(rep, cur_net, i + 1, self.tag, mtype.name + "_M")
                gmm.measure = mz
                gmm.found = count
                self.session.add(gmm)

                gmb = GeneMeasureRecord(rep, cur_net, i + 1, self.tag, mtype.name + "_B")
                gmb.measure = mz * rz
                gmb.found = count
                self.session.add(gmb)

                gmc = GeneMeasureRecord(rep, cur_net, i + 1, self.tag, mtype.name + "_C")
                gmc.measure = mz if mz < rz else rz
                gmc.found = count
                self.session.add(gmc)


    def analyse_multiple(self, rep, nets):
        for cur_net in nets:
            same, nov, dead = get_novel(
                cur_net, self.explore_samples, self.valid_patterns, self.optimal_rates)
            same_count = count_genes(cur_net, same)
            nov_count = count_genes(cur_net, nov)
            dead_count = count_genes(cur_net, dead)
            log.info("{} -- {}".format(cur_net, nov_count))

            for i in range(self.regs):
                rz = self.rz_info[i]
                mz = self.mz_info[i]
                for count, kind in zip((same_count, nov_count, dead_count), 
                                       "SAME NOVEL DEAD".split()):
                    gene_counts, total = count
                    if i in gene_counts:
                        gm = GeneMeasureRecord(rep, cur_net, i + 1, self.tag, kind)
                        gm.measure = rz
                        gm.found = gene_counts[i]
                        self.session.add(gm)

                        # gm.found = gene_counts[i] / float(total)
                        # gm.measure = info[i]



