import logging
log = logging.getLogger("core")

from enum import IntEnum
from .core_ext import SelectionModel, World, Channels, Population
from .analysis_ext import NetworkAnalysis
from .targets_ext import DefaultTarget, NoisyTarget

__all__ = ["World", "SelectionModel", "Parameters", "DefaultTarget", "NoisyTarget",
           "Channels", "NetworkAnalysis", "Population"]


class InterventionState(IntEnum):
    INTERVENE_NONE = 0
    INTERVENE_ON = 1
    INTERVENE_OFF = 2


class MutateType(IntEnum):
    JUMP = 0
    PROGRESSIVE = 1


class ScoringMethod(IntEnum):
    LINEAR = 0
    EXPONENTIAL = 1
    EXPONENTIAL_VEC = 2


class InputType(IntEnum):
    CONSTANT = 0
    PULSE = 1


class Parameters(object):

    def __init__(self, **kwargs):
        # Defaults are provided here
        self.seed = 1
        self.cis_count = 3
        self.gene_count = 3
        self.cue_channels = 2
        self.reg_channels = 1
        self.out_channels = 1
        self.selection_class = SelectionModel
        self.population_size = 100
        self.mutation_rate = 0.001
        self.trans_mutation_rate = 0.0
        self.duplication_rate = 0.0
        self.mutate_type = MutateType.JUMP
        self.score_method = ScoringMethod.LINEAR
        self.score_strength = 1.0
        self.add_zeros = 0
        self.input_type = InputType.CONSTANT
        self.pulse_for = 1

        self._override(kwargs)

    def _override(self, kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

        if not hasattr(self, 'reg_gene_count'):
            self.reg_gene_count = self.reg_channels

    def same_as(self, other):
        sdict = self.__dict__.copy()
        odict = other.__dict__.copy()

        # For now!
        for d in sdict, odict:
            if 'score_method' not in d:
                d['score_method'] = ScoringMethod.LINEAR
            if 'score_strength' not in d:
                d['score_strength'] = 1.0
            if 'input_type' not in d:
                d['input_type'] = InputType.CONSTANT
            if 'pulse_for' not in d:
                d['pulse_for'] = 1

        if set(sdict.keys()) != set(odict.keys()):
            x = set(sdict.keys())
            y = set(odict.keys())
            log.warning(
                "Differences in parameters-- {}/{}".format(x - y, y - x))
            return False

        for k, v in sdict.items():
            if odict[k] != v:
                log.warning("parameters {}, has changed: {} -> {}".format(
                    k, v, odict[k]))
                return False

        return True
