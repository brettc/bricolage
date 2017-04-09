import analysis_ext
import core_ext
import numpy

class NetworkAnalysis(analysis_ext.NetworkAnalysis):
    def __init__(self, network):
        super(NetworkAnalysis, self).__init__(network)
        self.annotations = {}

    def annotate(self, target):
        if target is not None:
            self.calc_mutual_info(target)
            self.calc_output_control(target)

    def calc_mutual_info(self, target):
        """Annotate the networks with information"""
        cats = target.calc_categories()
        if len(set(cats)) == 1:
            return

        analyzer = analysis_ext.MutualInfoAnalyzer(
            self.network.factory.world, cats)

        info = analyzer.numpy_info_from_network(self.network)
        info.shape = info.shape[1],
        self.mutual_info = info
        for i, mi in enumerate(info):
            a = self.annotations.setdefault(i, {})
            a['M'] = mi

    def calc_output_control(self, target):
        analyzer = analysis_ext.RelevantControlAnalyzer(
            self.network.factory.world,
            target.calc_distinct_outputs(),
            use_natural=False)
        info = analyzer.numpy_info_from_network(self.network)
        self.relevant_control = info

        for i, ri in enumerate(info):
            a = self.annotations.setdefault(i, {})
            # Make it available via the channel too
            a['R'] = ri


# WIP
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
        assert isinstance(collection, core_ext.CollectionBase)
        self._result = result
        self._collection = collection
        self._world = self._collection.factory.world
        self._array = numpy.asarray(self._result)
        self._output_size = self._world.out_channels

    @property
    def control(self):
        return self._array[:, :, :self._output_size]

    @property
    def entropy(self):
        return self._array[:, :, self._output_size:]

    @property
    def with_control(self):
        controlled = []
        for i, (c_net, e_net) in enumerate(zip(self.control, self.entropy)):
            for j, (c_reg, e_reg) in enumerate(zip(c_net, e_net)):
                # It must have full information about the environment
                if e_reg.prod() != 0.0:
                    if (e_reg == c_reg).all():
                        controlled.append(i)
        return controlled

    # @property
    # def fast_with_control(self):
    #     non_zero = self.entropy.prod(axis=2)


class AverageControlAnalyzer(analysis_ext.AverageControlAnalyzer):
    def __init__(self, world):
        super(AverageControlAnalyzer, self).__init__(world)

    def calc_info(self, net_or_collection):
        if isinstance(net_or_collection, core_ext.Network):
            res = self.analyse_network(net_or_collection)
            return AverageControlNetwork(net_or_collection, res)

        res = self.analyse_collection(net_or_collection)
        return AverageControlCollection(net_or_collection, res)
