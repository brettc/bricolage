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
        self.labeller = None

    def get_label(self, node_type, ident):
        if self.labeller is not None:
            return self.labeller.get_label(node_type, ident)

        if node_type == NodeType.GENE:
            return self.graph.get_gene_label(ident)
        elif node_type == NodeType.MODULE:
            return self.graph.get_module_label(ident)
        elif node_type == NodeType.CHANNEL:
            return self.graph.get_channel_label(ident)

        raise RuntimeError

    def get_node_attributes(self, node):
        name = self.graph.node_to_name(node)
        ntype, ident = node
        if ntype == NodeType.GENE:
            shape = 'hexagon' if self.graph.is_structural(node) else 'box'
            attrs = {
                'shape': shape,
                'label': self.get_label(ntype, ident),
            }
        elif ntype == NodeType.MODULE:
            attrs = {
                'shape': 'oval',
                'label': self.get_label(ntype, ident),
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
                'label': self.get_label(ntype, ident),
            }
        else:
            attrs = {}

        return name, attrs

    def get_dot(self, labeller=None):
        self.labeller = labeller

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

    def __init__(self, graph, height=8.0, width=5.0):
        super(DotDiagram, self).__init__()

        assert isinstance(graph, BaseGraph)
        self.graph = graph
        if hasattr(graph.analysis, 'annotations'):
            self.annotations = graph.analysis.annotations
        else:
            self.annotations = None

        self.world = graph.analysis.world

        self.width = width
        self.height = height

        # Set these below
        self.xscaling = 1.0
        self.yscaling = 1.0

        self.connections = []
        self.shapes_by_name = {}

        self._generate()

    def _generate(self):
        """Use dot program to initialise the diagram
        """

        if hasattr(self, 'get_label'):
            labeller = self
        else:
            labeller = None
        dot = DotMaker(self.graph).get_dot(labeller)

        # Lay it out
        dot.layout(prog='dot', args=_dot_default_args)

        sx, sy = self._calc_size(dot)
        self.yscaling = self.height / sy
        if self.width is not None:
            self.xscaling = float(self.width) / sx
        else:
            self.xscaling = self.yscaling

        for anode in dot.nodes_iter():
            self.add_shape(anode)

        for e in dot.edges():
            self.add_connection(e)

    def _calc_size(self, dot):
        """Calculate the size from the layout of the nodes
        """
        minx = miny = maxx = maxy = None

        def _update(mn, mx, cur):
            if mn is None or cur < mn:
                mn = cur
            if mx is None or cur > mx:
                mx = cur
            return mn, mx

        for anode in dot.nodes_iter():
            px, py = self.get_pt(anode.attr['pos'])
            minx, maxx = _update(minx, maxx, px)
            miny, maxy = _update(miny, maxy, py)

        return maxx - minx, maxy - miny


    def add_shape(self, anode):
        node_id = self.graph.name_to_node(anode)
        ntype, ident = node_id
        px, py = self.make_pt(anode.attr['pos'])

        if ntype == NodeType.GENE:
            gene = self.graph.network.genes[ident]
            shape = self.get_gene_shape(px, py, gene)
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

    def get_pt(self, strpair):
        x, y = strpair.split(',')
        return float(x), float(y)

    def make_pt(self, strpair):
        x, y = self.get_pt(strpair)
        x *= self.xscaling
        y *= self.yscaling
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
        ana.annotate(target)
        g = get_graph_by_type(graph_type, ana, knockouts=simplify)
        d = DotMaker(g)
        if name is None:
            name = str(net.identifier)
        path = path / name
        d.save_picture(str(path.with_suffix('.png')))
        if with_dot:
            d.save_dot(str(path.with_suffix('.dot')))
        return path
