import pathlib
from pygraphviz import AGraph
from pyx_drawing import Diagram

from bricolage.graph_maker import (NodeType, BaseGraph, GraphType,
                                   _decode_module_id, get_graph_by_type)
from analysis import NetworkAnalysis

_dot_default_args = '-Nfontname="Helvetica-8"'


class DotMaker(object):
    """Convert an NXGraph into an AGraph."""

    def __init__(self, graph):
        """Used to make dot files from network using Graphviz
        """
        assert isinstance(graph, BaseGraph)
        self.graph = graph

    def get_node_attributes(self, node):
        name = self.graph.node_to_name(node)
        ntype, ident = node
        if ntype == NodeType.GENE:
            # color = 'green' if self.graph.is_inert(node) else 'black'
            shape = 'hexagon' if self.graph.is_structural(node) else 'box'
            attrs = {
                'shape': shape,
                'label': self.graph.get_gene_label(ident),
            }
        elif ntype == NodeType.MODULE:
            attrs = {
                'shape': 'oval',
                'label': self.graph.get_module_label(ident),
            }

        elif ntype == NodeType.CHANNEL:
            if self.graph.is_input(node):
                shape = 'invtriangle'
            elif self.graph.is_output(node):
                shape = 'triangle'
            else:
                shape = 'diamond'

            attrs = {
                'shape': shape,
                'label': self.graph.get_channel_label(ident),
            }
        else:
            attrs = {}

        return name, attrs

    def get_dot(self):
        a_graph = AGraph(directed=True)
        nx_graph = self.graph.nx_graph

        # TODO: Add some default stuff?
        # a_graph.graph_attr.update(N.graph.get('graph',{}))
        # a_graph.node_attr.update(N.graph.get('node',{}))
        # a_graph.edge_attr.update(N.graph.get('edge',{}))

        structural_nodes = []
        output_nodes = []
        input_nodes = []
        # First, add nodes
        for node in nx_graph.nodes():
            name, attrs = self.get_node_attributes(node)
            if self.graph.is_input(node):
                input_nodes.append(name)
            elif self.graph.is_structural(node):
                structural_nodes.append(name)
            elif self.graph.is_output(node):
                output_nodes.append(name)

            # Keep a reference to the original node
            a_graph.add_node(name, **attrs)

        # We need to add subgraphs to cluster stuff on rank
        sub = a_graph.add_subgraph(input_nodes, name='input')
        sub.graph_attr['rank'] = 'source'
        sub = a_graph.add_subgraph(structural_nodes, name='structural')
        sub.graph_attr['rank'] = 'same'
        sub = a_graph.add_subgraph(output_nodes, name='output')
        sub.graph_attr['rank'] = 'sink'

        # Now add edges
        for u, v, edgedata in nx_graph.edges_iter(data=True):
            attrs = {}
            a_graph.add_edge(self.graph.node_to_name(u),
                             self.graph.node_to_name(v),
                             **attrs)

        return a_graph

    def save_picture(self, f):
        a = self.get_dot()
        a.draw(f, prog='dot', args=_dot_default_args)

    def save_dot(self, f):
        a = self.get_dot()
        a.write(f)


class DotDiagram(Diagram):
    # Map the binding types to arrows

    def __init__(self, graph, xscale=1.0, yscale=1.0):
        super(DotDiagram, self).__init__()
        # All the drawing shapes
        self.graph = graph

        # These just come from fiddling
        self.xscaling = .020 * xscale
        self.yscaling = .020 * yscale

        # self.gene_shapes = {}
        # self.signal_shapes = {}
        self.connections = []
        self.shapes_by_name = {}

        self._generate()

    def _generate(self):
        """Use dot program to initialise the diagram
        """
        dot = DotMaker(self.graph).get_dot()

        # Lay it out
        dot.layout(prog='dot', args=_dot_default_args)

        for anode in dot.nodes_iter():
            self.add_shape(anode)

        for e in dot.edges():
            self.add_connection(e)

    def add_shape(self, anode):
        node_id = self.graph.name_to_node(anode)
        ntype, ident = node_id
        px, py = self.make_pt(anode.attr['pos'])

        if ntype == NodeType.GENE:
            gene = self.graph.network.genes[ident]
            shape = self.get_gene_shape(px, py,  gene)
            self.add_object(shape, zorder=1)

        elif ntype == NodeType.CHANNEL:
            channel_type = self.graph.get_channel_type(node_id)
            shape = self.get_signal_shape(px, py, ident, channel_type)
            self.add_object(shape, zorder=2)

        elif ntype == NodeType.MODULE:
            gi, mi = _decode_module_id(ident)
            mod = self.graph.network.genes[gi].modules[mi]
            shape = self.get_module_shape(px, py, mod)
            self.add_object(shape, zorder=2)

        else:
            raise RuntimeError("Node type not found")

        self.shapes_by_name[anode] = shape

    def add_connection(self, edge):
        """Decipher the information loaded into the edge by the dot program"""
        n1, n2 = edge

        # Drop the 'e,' from the description, and read in all of the
        # points into an array
        point_string = edge.attr['pos'][2:]
        all_points = [self.make_pt(p) for p in point_string.split()]

        # We begin at point 1, and then skip each 3. Basically we're just
        # looking for the points describing the arc, not the control
        # points, as we reconstruct these.
        current_index = 1
        len_points = len(all_points)
        points = []
        while current_index < len_points:
            points.append(all_points[current_index])
            current_index += 3

        # We now have enough to make a nice curvy path
        shape1 = self.shapes_by_name[n1]
        shape2 = self.shapes_by_name[n2]

        connector = self.get_connector(shape1, shape2, points)
        self.add_object(connector, zorder=3)

    def make_pt(self, strpair):
        x, y = strpair.split(',')
        x = float(x) * self.xscaling
        y = float(y) * self.yscaling
        return x, y

    def get_gene_shape(self, px, py, gene):
        raise NotImplementedError

    def get_signal_shape(self, px, py, channel, channel_type):
        raise NotImplementedError

    def get_module_shape(self, px, py, mod):
        raise NotImplementedError

    def get_connector(self, shape1, shape2, points):
        raise NotImplementedError


def save_network_as_fullgraph(net, path='.', name=None,
                              simplify=True,
                              graph_type=GraphType.GENE_SIGNAL,
                              target=None,
                              with_dot=False):
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)
        ana = NetworkAnalysis(net)
        ana.calc_output_control()
        if target:
            ana.calc_mutual_info(target)
        g = get_graph_by_type(graph_type, ana)
        d = DotMaker(g)
        if name is None:
            name = str(net.identifier)
        path = path / name
        d.save_picture(str(path.with_suffix('.png')))
        if with_dot:
            d.save_dot(str(path.with_suffix('.dot')))
        return path
