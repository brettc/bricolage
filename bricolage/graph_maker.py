import networkx as nx
from enum import IntEnum

from .cis_logic import text_for_gene, text_for_cis_mod


class NodeType(IntEnum):
    # NOTE: These should be the same as those defined in pubsub2_c.h
    GENE = 0
    MODULE = 1
    CHANNEL = 2
    BEGIN = 3
    END = 4


class GraphType(IntEnum):
    FULL = 0
    SIGNAL_FLOW = 1
    GENE_SIGNAL = 2
    GENE = 3


class ChannelType(IntEnum):
    INPUT = 0
    INTERNAL = 1
    OUTPUT = 2


# NOTE: Should mirror "encode_module_id" in pubsub2_c.h
def _decode_module_id(module_id):
    return (0xff00 & module_id) >> 8, 0xff & module_id


def _encode_module_id(gene_id, module_id):
    return gene_id << 8 | module_id


class BaseGraph(object):
    def __init__(self, analysis, knockouts=True):
        self.analysis = analysis
        self.knockouts = knockouts
        self.original_network = analysis.network
        self.knockout_network = analysis.modified
        self.world = self.network.factory.world
        self.nx_graph = nx.DiGraph()

    @property
    def network(self):
        return self.knockout_network if self.knockouts else \
            self.original_network

    def is_inert(self, node):
        return False

    def get_channel_type(self, node):
        t, n = node
        assert t == NodeType.CHANNEL
        if self.is_input(node):
            return ChannelType.INPUT
        if self.is_internal(node):
            return ChannelType.INTERNAL
        return ChannelType.OUTPUT

    def is_internal(self, node):
        return not self.is_input(node) and not self.is_output(node)

    def is_input(self, node):
        t, n = node
        if t == NodeType.CHANNEL and n in self.world.cue_signals:
            return True
        return False

    def is_structural(self, node):
        t, n = node
        if t == NodeType.GENE and n >= self.world.reg_channels:
            return True
        return False

    def is_output(self, node):
        t, n = node
        if t == NodeType.CHANNEL and n in self.world.out_signals:
            return True
        return False

    def get_gene_description(self, gene_index):
        return str(self.network.genes[gene_index])

    def node_to_name(self, node):
        """Convert the node descriptor into a readable name"""
        ntype, nindex = node
        if ntype == NodeType.GENE:
            return "G{}".format(nindex)
        elif ntype == NodeType.MODULE:
            gindex, mindex = _decode_module_id(nindex)
            return "M{}-{}".format(gindex, mindex)
        elif ntype == NodeType.CHANNEL:
            return "C{}".format(nindex)
        elif ntype == NodeType.BEGIN:
            return "begin-{}".format(nindex)
        elif ntype == NodeType.END:
            return "end-{}".format(nindex)
        raise RuntimeError("Unknown node type {}".format(ntype))

    _ntype_lookup = {
        "G": NodeType.GENE, "M": NodeType.MODULE, "C": NodeType.CHANNEL
    }

    def name_to_node(self, name):
        """Convert the name back into a node descriptor"""
        ntype = self._ntype_lookup[name[:1]]
        if ntype == NodeType.MODULE:
            gindex, mindex = map(int, name[1:].split("-"))
            nindex = _encode_module_id(gindex, mindex)
        else:
            nindex = int(name[1:])
        return ntype, nindex

    def remove_nodes(self, nodetype, internal_only=False, external_only=False):
        assert not (internal_only and external_only)
        gr = self.nx_graph
        for nd in gr.nodes():
            if nd[0] != nodetype:
                continue

            if internal_only:
                if not self.is_internal(nd):
                    continue
            elif external_only:
                if self.is_internal(nd):
                    continue

            # Ok -- do the removal
            pred = gr.predecessors(nd)
            succ = gr.successors(nd)
            for p in pred:
                for s in succ:
                    gr.add_edge(p, s)
            gr.remove_node(nd)

    def format_annotations(self, obj):
        try:
            a = self.analysis.annotations[obj]
        except (AttributeError, KeyError):
            return ""
        else:
            keys = a.keys()
            keys.sort()
            return ", ".join(
                ["{}:{}".format(k, round(a[k], 2))
                 for k in keys])


class FullGraph(BaseGraph):
    def __init__(self, analysis, knockouts=True):
        BaseGraph.__init__(self, analysis, knockouts)
        # Build the nx_graph
        if self.knockouts:
            edges = analysis.get_active_edges()
        else:
            edges = analysis.get_edges()

        for nfrom, nto in edges:
            self.nx_graph.add_edge(nfrom, nto)

    def get_gene_label(self, i):
        # mods = self.network.genes[i].modules
        return "G{}".format(i + 1)

    def get_module_label(self, i):
        gi, mi = _decode_module_id(i)
        m = self.network.genes[gi].modules[mi]
        return text_for_cis_mod(self.world, m)

    def get_channel_label(self, i):
        return self.world.name_for_channel(i)


class GeneSignalGraph(FullGraph):
    def __init__(self, analysis, knockouts=True):
        FullGraph.__init__(self, analysis, knockouts)
        self.remove_nodes(NodeType.MODULE)

    def get_gene_label(self, i):
        glabel = FullGraph.get_gene_label(self, i)
        gene = self.network.genes[i]
        equation = text_for_gene(self.world, gene)
        ann = self.format_annotations(i)
        return "{} : {} [{}]".format(glabel, equation, ann)


class GeneGraph(GeneSignalGraph):
    def __init__(self, analysis, knockouts=True):
        GeneSignalGraph.__init__(self, analysis, knockouts)
        self.remove_nodes(NodeType.CHANNEL, internal_only=True)

    def get_gene_label(self, i):
        glabel = FullGraph.get_gene_label(self, i)
        g = self.network.genes[i]
        equation = text_for_gene(self.world, g)
        w = self.network.factory.world
        return "{}: {} => {}".format(glabel, equation,
                                     w.name_for_channel(g.pub))


class SignalFlowGraph(FullGraph):
    begin_node = (NodeType.BEGIN, 0)
    end_node = (NodeType.END, 0)

    def __init__(self, analysis):
        FullGraph.__init__(self, analysis, knockouts=True)
        gr = self.nx_graph

        inp_nodes = [n for n in gr.nodes_iter() if self.is_input(n)]
        for n in inp_nodes:
            gr.add_edge(self.begin_node, n)

        out_nodes = [n for n in gr.nodes_iter() if self.is_output(n)]
        for n in out_nodes:
            gr.add_edge(n, self.end_node)

        # Now remove the other stuff
        self.remove_nodes(NodeType.MODULE)
        self.remove_nodes(NodeType.GENE)

        # We've replaced the outside channels with the begin/end nodes
        self.remove_nodes(NodeType.CHANNEL, external_only=True)

    def minimum_cut(self):
        # Sometimes begin nodes may not even be in the graph!
        if self.begin_node not in self.nx_graph.nodes():
            return None
        return nx.minimum_node_cut(
            self.nx_graph, self.begin_node, self.end_node)


_type_map = {
    GraphType.FULL: FullGraph,
    GraphType.GENE_SIGNAL: GeneSignalGraph,
    GraphType.GENE: GeneGraph,
    GraphType.SIGNAL_FLOW: SignalFlowGraph,
}


def get_graph_by_type(gtype, *args, **kwargs):
    return _type_map[gtype](*args, **kwargs)