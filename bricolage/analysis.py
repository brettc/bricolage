import analysis_ext
import core_ext
import numpy


class NetworkAnalysis(analysis_ext.NetworkAnalysis):
    def __init__(self, network):
        super(NetworkAnalysis, self).__init__(network)
        self.annotations = {}

    def calc_mutual_info(self, target):
        """Annotate the networks with information"""
        analyzer = analysis_ext.MutualInfoAnalyzer(
            self.network.factory.world,
            target.calc_categories())

        # TODO: This should be done in the function!
        info = analyzer.numpy_info_from_network(self.network)
        info.shape = info.shape[1],
        self.mutual_info = info

        for mi, g in zip(info, self.network.genes):
            a = self.annotations.setdefault(g.sequence, {})
            a['M'] = mi

    def calc_output_control(self):
        analyzer = analysis_ext.OutputControlAnalyzer(self.network.factory.world)
        info = analyzer.numpy_info_from_network(self.network)
        info.shape = info.shape[1:]
        self.output_control = info

        for oi, g in zip(info, self.network.genes):
            a = self.annotations.setdefault(g.sequence, {})
            # Make it available via the channel too
            a['C'] = oi[0]
            a['E'] = oi[1]


# TODO: Possible way of handling information in more pythonic way...
class AverageControlNetwork(object):
    def __init__(self, net, result):
        self._result = result
        self._net = net
        arr = numpy.asarray(result)
        assert arr.shape[0] == 1
        arr.shape = arr.shape[1:]
        self._array = arr


class AverageControlCollection(object):
    def __init__(self, collection, result):
        self._result = result
        self._collection = collection


class AverageControlAnalyzer(analysis_ext.AverageControlAnalyzer):
    def __init__(self, world):
        super(AverageControlAnalyzer, self).__init__(world)

    def calc_info(self, net_or_collection):
        if isinstance(net_or_collection, core_ext.Network):
            res = self.analyse_network(net_or_collection)
            return AverageControlNetwork(net_or_collection, res)

        res = self.analyse_collection(net_or_collection)
        return AverageControlCollection(net_or_collection, res)
