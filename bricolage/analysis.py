import analysis_ext


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
