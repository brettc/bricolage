"""Useful functions for loading pandasa frames into ipython notebook
"""

import numpy as np
import pandas as pd

from .analysis_ext import (
    MutualInfoAnalyzer, CausalFlowAnalyzer, AverageControlAnalyzer,
    Information, NetworkAnalysis)
from .graph import SignalFlowGraph
from neighbourhood import NetworkNeighbourhood


def make_population_frames(pop, target, flow, do_cuts=True):
    # Assumption -- fitness is calculated!
    # Create the fitnesses
    fits = np.zeros(pop.size)
    pop.get_fitnesses(fits)
    fits_frame = pd.DataFrame({'fitness': fits})

    # The causal flow analysis
    cf = CausalFlowAnalyzer(pop.factory.world, flow)
    ci = cf.numpy_info_from_collection(pop)
    csummed = ci.sum(axis=2)
    cframe = pd.DataFrame(csummed)
    cframe.columns = ["C{}".format(i) for i in range(1, len(cframe.columns) + 1)]

    # Control analysis
    af = AverageControlAnalyzer(pop.factory.world, flow)
    ai = np.asarray(af.analyse_collection(pop))
    asummed = ai.sum(axis=2)
    aframe = pd.DataFrame(asummed)
    aframe.columns = ["A{}".format(i) for i in range(1, len(aframe.columns) + 1)]

    mf = MutualInfoAnalyzer(pop.factory.world, target.calc_categories())
    mi = mf.numpy_info_from_collection(pop)
    msummed = mi.sum(axis=2)
    mframe = pd.DataFrame(msummed)
    mframe.columns = ["M{}".format(i) for i in range(1, len(mframe.columns) + 1)]

    # Cuts are time-consuming...
    if not do_cuts:
        return fits_frame, mframe, cframe, aframe

    cuts = []
    for n in pop:
        ana = NetworkAnalysis(n)
        fg = SignalFlowGraph(ana)
        cuts.append(len(fg.minimum_cut()))
    cuts_frame = pd.DataFrame({'cuts': cuts})

    return fits_frame, cuts_frame, mframe, cframe, aframe


def make_population_frame(pop, target, flow, do_cuts=True):
    return pd.concat(make_population_frames(pop, target, flow, do_cuts), axis=1)


def make_ancestry_frames(anc, target, flow):
    target.assess_collection(anc)

    # Create a generations dataframe for indexing
    # gens = pd.DataFrame({'generation': [n.generation for n in anc]})
    gens = np.asarray([n.generation for n in anc], dtype=np.int64)

    # Create the fitnesses
    fits_frame = pd.DataFrame({'fitness': [n.fitness for n in anc]})
    fits_frame.index = gens

    # The causal flow analysis
    cf = CausalFlowAnalyzer(anc.factory.world, flow)
    cj = cf.analyse_collection(anc)
    ci = np.asarray(Information(cj))
    csummed = ci.sum(axis=2)
    cframe = pd.DataFrame(csummed)
    cframe.index = gens
    cframe.columns = ["C{}".format(i) for i in range(1, len(cframe.columns) + 1)]

    # The mutual flow analysis
    mf = MutualInfoAnalyzer(anc.factory.world, target.calc_categories())
    mj = mf.analyse_collection(anc)
    mi = np.asarray(Information(mj))
    msummed = mi.sum(axis=2)
    mframe = pd.DataFrame(msummed)
    mframe.index = gens
    mframe.columns = ["M{}".format(i) for i in range(1, len(mframe.columns) + 1)]

    return fits_frame, mframe, cframe


def get_population_neighbourhood_fitness(pop, target, sample_per_network=1000):
    fits = np.zeros(sample_per_network)

    means = []

    for n in pop:
        nayb = NetworkNeighbourhood(n, sample_per_network)
        target.assess_collection(nayb.neighbours)
        nayb.neighbours.get_fitnesses(fits)
        means.append(fits.mean())

    npmeans = np.array(means)
    return npmeans.mean()


def make_ancestry_robustness_frame(anc, target, sample_size=1000):
    # Get the fits of the ancestry
    target.assess_collection(anc)
    fits = np.zeros(anc.size)
    anc.get_fitnesses(fits)

    # Now sample the region around each individual and get the mean fitness
    sample_fits = np.zeros(sample_size)
    means = np.zeros(anc.size)
    top_count = np.zeros(anc.size)
    for i, n in enumerate(anc):
        nayb = NetworkNeighbourhood(n, sample_size)
        target.assess_collection(nayb.neighbours)
        nayb.neighbours.get_fitnesses(sample_fits)
        how_many_are_1 = (sample_fits == 1.0).sum()
        means[i] = sample_fits.mean()
        top_count[i] = float(how_many_are_1) / sample_size

    df = pd.DataFrame({'fitness': fits, 'n_mean': means, 'max_count': top_count})
    gens = np.asarray([n.generation for n in anc], dtype=np.int64)
    df.index = gens

    return df


def make_cut_frame(anc):
    gens = []
    cuts = []

    for n in anc:
        gens.append(n.generation)
        ana = NetworkAnalysis(n)
        fg = SignalFlowGraph(ana)
        cuts.append(len(fg.minimum_cut()))

    df = pd.DataFrame({'cuts': cuts})
    df.index = gens
    return df
