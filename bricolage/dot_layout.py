import pathlib
from pygraphviz import AGraph
from pyx_drawing import Diagram
from bricolage.graph_maker import (NodeType, BaseGraph, GraphType,
                                   decode_module_id, get_graph_by_type,
                                   node_logic_differs)
from analysis import NetworkAnalysis

_dot_default_args = '-Nfontname="Helvetica-8"'


class DotMaker(object):
    """Convert an NXGraph into an AGraph."""

    def __init__(self, graph, simple=False):
        """Used to make dot files from network using Graphviz
        """
        assert isinstance(graph, BaseGraph)
        self.graph = graph
        self.labeller = None
        self.simple = simple

    # def get_label(self, graph, node_type, ident):
    #     if self.labeller is not None:
    #         return self.labeller.get_label(node_type, ident)
    #
    #     return graph.get_label((node_type, ident))

    def get_node_attributes(self, node, use_graph=None):
        if use_graph is None:
            use_graph = self.graph

        name = use_graph.node_to_name(node)
        ntype, ident = node
        if ntype == NodeType.GENE:
            shape = 'hexagon' if use_graph.is_structural(node) else 'box'
            attrs = {
                'shape': shape,
                'label': use_graph.get_label((ntype, ident), self.simple),
            }
        elif ntype == NodeType.MODULE:
            attrs = {
                'shape': 'oval',
                'label': use_graph.get_label((ntype, ident), self.simple),
            }

        elif ntype == NodeType.CHANNEL:
            if use_graph.is_input(node):
                shape = 'oval'
                # shape = 'invtriangle'
            elif use_graph.is_output(node):
                shape = 'triangle'
            else:
                shape = 'diamond'

            attrs = {
                'shape': shape,
                'label': use_graph.get_label((ntype, ident), self.simple),
            }
        else:
            attrs = {}

        return name, attrs

    def categorize_node(self, node, name,
                        input_nodes, structural_nodes, output_nodes):
        if self.graph.is_input(node):
            input_nodes.append(name)
        elif self.graph.is_structural(node):
            structural_nodes.append(name)
        elif self.graph.is_output(node):
            output_nodes.append(name)

    def get_dot(self, labeller=None):
        self.labeller = labeller

        a_graph = AGraph(directed=True)
        nx_graph = self.graph.nx_graph

        structural_nodes = []
        output_nodes = []
        input_nodes = []
        # First, add nodes
        for node in nx_graph.nodes():
            name, attrs = self.get_node_attributes(node)

            self.categorize_node(node, name, input_nodes, structural_nodes, output_nodes)

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
        for u, v, edgedata in nx_graph.edges(data=True):
            attrs = {}
            a_graph.add_edge(self.graph.node_to_name(u),
                             self.graph.node_to_name(v),
                             **attrs)

        return a_graph

    def get_dot_diff(self, other):

        a_graph = AGraph(directed=True)
        nx_to_graph = self.graph.nx_graph
        nx_from_graph = other.nx_graph

        structural_nodes = []
        output_nodes = []
        input_nodes = []

        # First, Add the nodes from the "to" Graph, marking those that are
        # changed and new.
        for node in nx_to_graph.nodes():
            name, attrs = self.get_node_attributes(node)

            if node in nx_from_graph.nodes():
                if node_logic_differs(self.graph, other, node):
                    attrs['color'] = 'blue'
            else:
                attrs['color'] = 'green'

            self.categorize_node(node, name,
                                 input_nodes, structural_nodes, output_nodes)

            # Keep a reference to the original node
            a_graph.add_node(name, **attrs)

        for node in nx_from_graph.nodes():
            name, attrs = self.get_node_attributes(node, use_graph=other)
            if node not in nx_to_graph.nodes():
                attrs['color'] = 'red'
                self.categorize_node(node, name,
                                     input_nodes, structural_nodes, output_nodes)
                a_graph.add_node(name, **attrs)


        # We need to add subgraphs to cluster stuff on rank
        sub = a_graph.add_subgraph(input_nodes, name='input')
        sub.graph_attr['rank'] = 'source'
        sub = a_graph.add_subgraph(structural_nodes, name='structural')
        sub.graph_attr['rank'] = 'same'
        sub = a_graph.add_subgraph(output_nodes, name='output')
        sub.graph_attr['rank'] = 'sink'

        # Now add edges
        for u, v, edgedata in nx_to_graph.edges(data=True):
            attrs = {}
            if (u, v) not in nx_from_graph.edges():
                attrs['color'] = 'green'

            a_graph.add_edge(self.graph.node_to_name(u),
                             self.graph.node_to_name(v),
                             **attrs)

        for edge in nx_from_graph.edges():
            if edge not in nx_to_graph.edges():
                attrs = {'color': 'red'}
                a_graph.add_edge(self.graph.node_to_name(edge[0]),
                                self.graph.node_to_name(edge[1]),
                                **attrs)


        return a_graph

    def save_picture(self, f):
        a = self.get_dot()
        a.draw(f, prog='dot', args=_dot_default_args)

    def save_dot(self, f):
        a = self.get_dot()
        a.write(f)

    def save_diff_picture(self, f, other):
        a = self.get_dot_diff(other)
        a.draw(f, prog='dot', args=_dot_default_args)

    def save_diff_dot(self, f, other):
        a = self.get_dot_diff(other)
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

        for anode in dot.nodes():
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

        for anode in dot.nodes():
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
            gi, mi = decode_module_id(ident)
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


def save_graph(net, path='.', name=None,
                 simplify=True,
                 graph_type=GraphType.GENE_SIGNAL,
                 target=None,
                 with_dot=False,
                 diff=None,
                 simple_labels=False):
    """Put it all together into a simple call
    """
    if not isinstance(path, pathlib.Path):
        path = pathlib.Path(path)
    ana = NetworkAnalysis(net)
    ana.annotate(target)
    g = get_graph_by_type(graph_type, ana, knockouts=simplify)
    d = DotMaker(g, simple=simple_labels)
    if name is None:
        name = str(net.identifier)
    path = path / name
    if diff is None:
        d.save_picture(str(path.with_suffix('.png')))
        if with_dot:
            d.save_dot(str(path.with_suffix('.dot')))
    else:
        other_ana = NetworkAnalysis(diff)
        other_ana.annotate(target)
        gdiff = get_graph_by_type(graph_type, other_ana, knockouts=simplify)
        d.save_diff_picture(str(path.with_suffix('.png')), other=gdiff)
        if with_dot:
            d.save_diff_dot(str(path.with_suffix('.dot')), other=gdiff)

    return path
