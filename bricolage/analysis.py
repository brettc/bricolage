"""High level analysis tools
"""

import numpy as np
import pandas as pd

from analysis_ext import MutualInfoAnalyzer, CausalFlowAnalyzer, Information
from lineage import FullLineage

class LineageSummarizer(object):
    def __init__(self, lin, target, flow):
        assert isinstance(lin, FullLineage)
        self._lineage = lin
        self._cf = CausalFlowAnalyzer(lin.world, flow)
        self._mf = MutualInfoAnalyzer(lin.world, target.calc_categories())
        self.params = lin.params
        self.target = target
        self.size = lin.generation + 1
        self.fitnesses = np.zeros(self.params.population_size)

        self.data = np.zeros(self.size, self.dtype())

        # Allocate the stuff
    def dtype(self):
        regs = self.params.reg_channels

        desc = [
            ('F_25', np.double),
            ('F_50', np.double),
            ('F_75', np.double),

            ('C_25', (np.double, regs)),
            ('C_50', (np.double, regs)),
            ('C_75', (np.double, regs)),
        ]

        return np.dtype(desc)

    def calc_generations(self):
        for i, g in self._lineage.all_generations():
            self._calc(i, g)
            yield i, g

    def _calc(self, i, g):
        fits = self.fitnesses
        g.get_fitnesses(fits)
        self.data['F_25'][i] = np.percentile(fits, 25)
        self.data['F_50'][i] = np.percentile(fits, 50)
        self.data['F_75'][i] = np.percentile(fits, 75)

        cj = self._cf.analyse_collection(g)
        ci = np.asarray(Information(cj))
        csummed = ci.sum(axis=2)
        self.data['C_25'][i] = np.percentile(csummed, 25)
        self.data['C_50'][i] = np.percentile(csummed, 50)
        self.data['C_75'][i] = np.percentile(csummed, 75)

    def get_frame(self):
        # Spread the frames out...
        d = {}
        for i in 25, 50, 75:
            s = 'F_' + str(i)
            d[s] = self.data[s]

        regs = self.params.reg_channels
        for c in range(regs):
            for i in 25, 50, 75:
                fr = 'C_{}'.format(i)
                to = 'C{}_{}'.format(c + 1, i)
                d[to] = self.data[fr][:, c]
        return pd.DataFrame(d)


def make_population_frames(pop, target, flow):
    # Assumption -- fitness is calculated!
    # Create the fitnesses
    fits_frame = pd.DataFrame({'fitness': [n.fitness for n in pop]})

    # The causal flow analysis
    cf = CausalFlowAnalyzer(pop.constructor.world, flow)
    cj = cf.analyse_collection(pop)
    ci = np.asarray(Information(cj))
    csummed = ci.sum(axis=2)
    cframe = pd.DataFrame(csummed)
    cframe.columns = ["C{}".format(i) for i in range(1, len(cframe.columns) + 1)]

    mf = MutualInfoAnalyzer(pop.constructor.world, target.calc_categories())
    mj = mf.analyse_collection(pop)
    mi = np.asarray(Information(mj))
    msummed = mi.sum(axis=2)
    mframe = pd.DataFrame(msummed)
    mframe.columns = ["M{}".format(i) for i in range(1, len(mframe.columns) + 1)]

    return fits_frame, mframe, cframe


def make_ancestry_frames(anc, target, flow):
    target.assess_collection(anc)

    # Create a generations dataframe for indexing
    gens = pd.DataFrame({'generation': [n.generation for n in anc]})

    # Create the fitnesses
    fits_frame = pd.DataFrame({'fitness': [n.fitness for n in anc]})
    fits_frame.index = gens

    # The causal flow analysis
    cf = CausalFlowAnalyzer(anc.constructor.world, flow)
    cj = cf.analyse_collection(anc)
    ci = np.asarray(Information(cj))
    csummed = ci.sum(axis=2)
    cframe = pd.DataFrame(csummed)
    cframe.index = gens
    cframe.columns = ["C{}".format(i) for i in range(1, len(cframe.columns) + 1)]

    # The mutual flow analysis
    mf = MutualInfoAnalyzer(anc.constructor.world, target.calc_categories())
    mj = mf.analyse_collection(anc)
    mi = np.asarray(Information(mj))
    msummed = mi.sum(axis=2)
    mframe = pd.DataFrame(msummed)
    mframe.index = gens
    mframe.columns = ["M{}".format(i) for i in range(1, len(mframe.columns) + 1)]

    return fits_frame, mframe, cframe
