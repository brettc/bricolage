import pathlib
from pygraphviz import AGraph

from bricolage.graph_maker import NodeType, FullGraph, BaseGraph


class DotMaker(object):
    """Convert an NXGraph into an AGraph."""

    _dot_default_args = '-Nfontname="Helvetica-8"'

    def __init__(self, graph):
        """Used to make dot files from network using Graphviz
        """
        assert isinstance(graph, BaseGraph)
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

    def get_dot(self):
        a_graph = AGraph(directed=True)
        nx_graph = self.graph.nx_graph

        # TODO: Add some default stuff?
        # a_graph.graph_attr.update(N.graph.get('graph',{}))
        # a_graph.node_attr.update(N.graph.get('node',{}))
        # a_graph.edge_attr.update(N.graph.get('edge',{}))
        #
        # TODO: We could add some nice clustering stuff here, by annotating
        # the graph with clustering or "modules".

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
            # attrs['gnode'] = node
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
        a.draw(f, prog='dot', args=self._dot_default_args)

    def save_dot(self, f):
        a = self.get_dot()
        a.write(f)

    def make_diagram(self, layout):
        """Use dot program to initialise the layout
        """
        dot = self.get_dot()

        # Now lay it out
        dot.layout(prog='dot', args=self._dot_default_args)

        # New we extract the information from the dot file to give us
        # coordinates
        for anode in dot.nodes_iter():
            px, py = layout.make_pt(anode.attr['pos'])
            node_id = self.graph.name_to_node(anode)
            data = anode.attr
            layout.add_shape(node_id, px, py, data)

        for e in dot.edges():
            n1, n2 = e
            data = self.aedge_to_graph_data(e)
            # Drop the 'e,' from the description, and read in all of the
            # points into an array
            point_string = e.attr['pos'][2:]
            all_points = [layout.make_pt(p) for p in point_string.split()]

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
            shape1 = self.aname_to_shape(layout, n1)
            shape2 = self.aname_to_shape(layout, n2)

            layout.add_connection(shape1, shape2, points, data)

    def aedge_to_graph_data(self, aedge):
        name1, name2 = aedge
        node1, node2 = map(self.graph.name_to_node, [name1, name2])
        data = self.graph.nx_graph.get_edge_data(node1, node2)
        return data

    def aname_to_shape(self, layout, aname):
        t, i = self.graph.name_to_node(aname)
        if t == NodeType.GENE:
            return layout.gene_shapes[i]
        elif t == NodeType.CHANNEL:
            return layout.signal_shapes[i]
        raise RuntimeError("Node name not found")


class BaseLayout(object):
    # Map the binding types to arrows

    def __init__(self, diagram):
        # All the drawing shapes
        self.diagram = diagram

        # These just come from fiddling
        self.xscaling = .020
        self.yscaling = .020

        self.gene_shapes = {}
        self.signal_shapes = {}
        self.connections = []

    def get_gene_description(self, g):
        # return "G%02d: %s" % (g.index, g.description)
        # return "%s" % (g.description)
        return "XXXXXXXX"

    def add_shape(self, node_id, px, py, data):
        # Use the same questions we do above.
        t, i = node_id
        if t == NodeType.GENE:
            s = self.get_gene_shape(node_id, px, py, data)
            self.gene_shapes[i] = s
        elif t == NodeType.CHANNEL:
            s = self.get_signal_shape(node_id, px, py, data)
            self.signal_shapes[i] = s
        else:
            raise RuntimeError("Node type not found")

    def add_connection(self, shape1, shape2, points, data):
        # TODO: Just defaulting to arrows for now
        c = self.get_connection(shape1, shape2, points, data)
        self.connections.append(c)

    def make_pt(self, strpair):
        x, y = strpair.split(',')
        x = float(x) * self.xscaling
        y = float(y) * self.yscaling
        return x, y

    def get_gene_shape(self, gene_number, px, py, data):
        raise NotImplementedError

    def get_signal_shape(self, signal_number, px, py, data):
        raise NotImplementedError

    def get_connection(self, shape1, shape2, points, data):
        raise NotImplementedError

    def draw(self, C):
        for do in self.gene_shapes.values():
            do.draw(C)
        for do in self.signal_shapes.values():
            do.draw(C)
        for do in self.connections:
            do.draw(C)


def save_network_as_fullgraph(n, path='.', name=None, simplify=True):
    output_path = pathlib.Path(path)
    ana = NetworkAnalysis(n)
    gph = FullGraph(ana, simplify)
    dot = DotMaker(gph)
    if name is None:
        name = n.identifier
    print 'saving', name
    dot.save_picture(str(output_path / "network-{}.png".format(name)))
    dot.save_dot(str(output_path / "network-{}.dot".format(name)))
