from enum import IntEnum
import networkx as nx
from pygraphviz import AGraph
from core_ext import NetworkAnalysis
import pathlib

class NodeType(IntEnum):
    # NOTE: These should be the same as those defined in pubsub2_c.h
    GENE = 0
    MODULE = 1
    CHANNEL = 2
    BEGIN = 3
    END = 4


# NOTE: Should mirror "encode_module_id" in pubsub2_c.h
def _decode_module_id(module_id):
    return (0xff00 & module_id) >> 8, 0xff & module_id


def _encode_module_id(gene_id, module_id):
    return gene_id << 8 | module_id


class BaseGraph(object):
    def __init__(self, analysis):
        self.analysis = analysis
        self.network = analysis.network
        self.world = self.network.constructor.world
        self.nx_graph = nx.DiGraph()

    def is_inert(self, node):
        return False

    # def is_internal(self, node):
    #     t, n = node
    #     if t == "G":
    #         if n < self.params.reg_gene_count:
    #             return True
    #     else:
    #         if n in self.params.reg_signals:
    #             return True
    #     return False

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

    _ntype_lookup = {"G":NodeType.GENE, "M":NodeType.MODULE, "C":NodeType.CHANNEL}
    def name_to_node(self, name):
        """Convert the name back into a node descriptor"""
        ntype = self._ntype_lookup[name[:1]]
        if ntype == NodeType.MODULE:
            gindex, mindex = map(int, name[1:].split("-"))
            nindex = _encode_module_id(gindex, mindex)
        else:
            nindex = int(name[1:])
        return ntype, nindex

    def remove_nodes(self, nodetype, internal_only=False):
        G = self.nx_graph
        for nd in G.nodes():
            if nd[0] != nodetype:
                continue
            if internal_only:
                if self.is_input(nd) or self.is_output(nd):
                    continue
            pred = G.predecessors(nd)
            succ = G.successors(nd)
            for p in pred:
                for s in succ:
                    G.add_edge(p, s)
            G.remove_node(nd)


class FullGraph(BaseGraph):
    def __init__(self, analysis, knockouts=True):
        BaseGraph.__init__(self, analysis)
        # Build the nx_graph
        if knockouts:
            edges = analysis.get_active_edges()
        else:
            edges = analysis.get_edges()

        for nfrom, nto in edges:
            self.nx_graph.add_edge(nfrom, nto)

    def get_gene_label(self, i):
        # mods = self.network.genes[i].modules
        # return "|".join([str(j) for j, m in enumerate(mods)])
        return str(self.network.genes[i])
    
    def get_module_label(self, i):
        gi, mi = _decode_module_id(i)
        m = self.network.genes[gi].modules[mi]
        return str(m)

    def get_channel_label(self, i):
        return self.world.name_for_channel(i)


class SignalFlowGraph(FullGraph):
    def __init__(self, analysis):
        FullGraph.__init__(self, analysis, knockouts=True)
        G = self.nx_graph
        self.begin_node = (NodeType.BEGIN, 0)
        self.end_node = (NodeType.END, 0)

        # Replace input nodes with "begin"
        inp_nodes = [n for n in G.nodes_iter() if self.is_input(n)]
        for n in inp_nodes:
            G.add_edge((NodeType.BEGIN, 0), n)
            succ = G.successors(n)
            for s in succ:
                G.add_edge(self.begin_node, s)
            G.remove_node(n)

        # Replace output nodes with "end"
        out_nodes = [n for n in G.nodes_iter() if self.is_output(n)]
        for n in out_nodes:
            pred = G.predecessors(n)
            for p in pred:
                G.add_edge(p, self.end_node)
            G.remove_node(n)

        # Now remove the other stuff
        self.remove_nodes(NodeType.MODULE)
        self.remove_nodes(NodeType.GENE)

    def minimum_cut(self):
        return nx.minimum_node_cut(
            self.nx_graph, self.begin_node, self.end_node)

class GeneSignalGraph(FullGraph):
    def __init__(self, analysis, knockouts=True):
        FullGraph.__init__(self, analysis, knockouts)
        self.remove_nodes(NodeType.MODULE)

class GeneGraph(FullGraph):
    def __init__(self, analysis, knockouts=True):
        FullGraph.__init__(self, analysis, knockouts)
        self.remove_nodes(NodeType.MODULE)
        self.remove_nodes(NodeType.CHANNEL, internal_only=True)

class DotMaker(object):
    """Convert an NXGraph into an AGraph."""
    # Map the binding types to arrows.
    # arrow_types = {
    #     BindingType.INERT: "dot",
    #     BindingType.ACTIVE: "normal",
    #     BindingType.REPRESS: "tee",
    #     BindingType.COMPLEX: "empty",
    # }
    _dot_default_args = '-Nfontname="Helvetica-8"'

    def __init__(self, graph):
        """Used to make dot files from network using Graphviz
        """
        self.graph = graph

    def get_node_attributes(self, node):
        name = self.graph.node_to_name(node)
        t, i = node
        if t == NodeType.GENE:
            # color = 'green' if self.graph.is_inert(node) else 'black'
            shape = 'hexagon' if self.graph.is_structural(node) else 'box'
            attrs = {
                'shape': shape,
                'label': self.graph.get_gene_label(i),
            }
        elif t == NodeType.MODULE:
            attrs = {
                'shape': 'oval',
                'label': self.graph.get_module_label(i),
            }

        elif t == NodeType.CHANNEL:
            if self.graph.is_input(node):
                shape = 'invtriangle'
            elif self.graph.is_output(node):
                shape = 'triangle'
            else:
                shape = 'diamond'

            attrs = {
                'shape': shape,
                'label': self.graph.get_channel_label(i),
            }
        else:
            attrs = {}

        return name, attrs

    def get_layout(self):
        A = AGraph(directed=True)
        G = self.graph.nx_graph

        # TODO: Add some default stuff?
        # A.graph_attr.update(N.graph.get('graph',{}))
        # A.node_attr.update(N.graph.get('node',{}))
        # A.edge_attr.update(N.graph.get('edge',{}))
        #
        # TODO: We could add some nice clustering stuff here, by annotating
        # the graph with clustering or "modules".

        # internal = [] 
        structural_nodes = [] 
        output_nodes = []
        input_nodes = []
        # First, add nodes
        for node in G.nodes():
            name, attrs = self.get_node_attributes(node)
            if self.graph.is_input(node):
                input_nodes.append(name)
            elif self.graph.is_structural(node):
                structural_nodes.append(name)
            elif self.graph.is_output(node):
                output_nodes.append(name)

            # if self.is_internal(node):
            #     internal.append(name)
            # if node[0] == "G" and node[1] >= self.params.reg_gene_count:
            
            # if self.graph.is_node_old(node):
            #     attrs['style'] = 'dotted'
            # elif self.graph.is_node_new(node):
            #     attrs['color'] = 'red'
            
            # Keep a reference to the original node
            # attrs['gnode'] = node
            A.add_node(name, **attrs)

        # NOTE: bug in graphviz, if you want clustering, you must have a name
        # that begins "cluster"
        # A.add_subgraph(internal, name='cluster1')
        sub = A.add_subgraph(input_nodes, name='input')
        sub.graph_attr['rank'] = 'source'
        sub = A.add_subgraph(structural_nodes, name='structural')
        sub.graph_attr['rank'] = 'same'
        sub = A.add_subgraph(output_nodes, name='output')
        sub.graph_attr['rank'] = 'sink'

        # Now add edges
        for u, v, edgedata in G.edges_iter(data=True):
            # k = edgedata['kind']
            # c = 'green' if k == BindingType.INERT else 'black'
            # a = self.arrow_types[k]
            # attrs = {
            #     'color': c,
            #     'arrowhead': a,
            # }
            # if self.graph.is_edge_old(u, v):
            #     attrs['style'] = 'dotted'
            # elif self.graph.is_edge_new(u, v):
            #     attrs['color'] = 'red'
            #
            attrs = {}
            A.add_edge(self.graph.node_to_name(u), 
                self.graph.node_to_name(v), 
                **attrs)

        return A

    def save_picture(self, f):
        a = self.get_layout()
        a.draw(f, prog='dot', args=self._dot_default_args)

    def save_dot(self, f):
        a = self.get_layout()
        a.write(f)

def save_network_as_fullgraph(n, path='.', name=None, simplify=True):
    output_path = pathlib.Path(path)
    ana = NetworkAnalysis(n)
    gph = FullGraph(ana, simplify)
    dot = DotMaker(gph)
    if name is None:
        name = n.identifier
    print 'saving', name
    dot.save_picture(str(output_path / "network-{}.png".format(name)))



