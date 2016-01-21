import logging

log = logging.getLogger("drawing")

from graph_maker import FullGraph, GeneSignalGraph
from graph_layout import DotMaker, BaseLayout
from pyx_drawing import (Diagram, rounded_rect, CurvePath, Shape, cut_hexagon,
                       hexagon, diamond)
from pyx import canvas, color as pcolor, path as ppath


class ColorScheme(object):
    sig = pcolor.cmyk.CarnationPink
    cue = pcolor.cmyk.Blue
    act = pcolor.cmyk.Green
    # gene = pcolor.cmyk.Black
    gene = pcolor.cmyk.Gray
    gene_old = pcolor.cmyk.Gray
    gene_new = pcolor.cmyk.Red


class GrayScheme(object):
    channel = pcolor.cmyk.Black
    gene = pcolor.cmyk.Black
    sig = pcolor.cmyk.White
    cue = pcolor.cmyk.White
    act = pcolor.cmyk.White


def signal_to_variable(p, c):
    if c in p.cue_signals:
        txt = "e$_{%d}$" % (c + 1)
    elif c in p.reg_signals:
        txt = "r$_{%d}$" % (c - p.reg_range[0] + 1)
    else:
        txt = "s$_{%d}$" % (c - p.out_range[0] + 1)
    return txt


def latexify(l):
    l = l.replace("~", r"$\lnot{}$")
    # l = l.replace("~", r"$\sim$")
    # l = l.replace("~", "NOT ")
    # l = l.replace("~", r"\~{}")
    l = l.replace(" AND ", r" $\land{}$ ")
    l = l.replace(" OR ", r" $\lor{}$ ")
    l = l.replace(")", "")
    l = l.replace("(", "")
    l = l.replace("|", r"$\to$")
    return l


class SimpleGeneShape(Shape):
    width = .8
    height = .2

    def __init__(self, diag, px, py, node_id, data):
        Shape.__init__(self, diag, px, py)
        self.outline = rounded_rect(px, py, self.width, self.height)
        self.node_id = node_id
        self.gene_number = node_id[1]
        # Convert from stupid attribute type
        self.data = dict([(k, v) for k, v in data.items()])
        # self.gene = self.diagram.net.genes[self.gene_number]

    def draw(self, C):
        # No text labeling, just filled in black
        # if self.diagram.graph.'old' in self.data:
        #     color = self.diagram.scheme.gene_old
        #     print 'old'
        # if 'new' in self.data:
        #     color = self.diagram.scheme.gene_new
        #     print 'new'
        # else:
        # color = self.diagram.scheme.gene
        C.fill(self.outline, [pcolor.cmyk.Gray])


class AGeneShape(Shape):
    width = .8
    height = .8

    def __init__(self, diag, px, py, node_id, data):
        Shape.__init__(self, diag, px, py)
        self.outline = rounded_rect(px, py, self.width, self.height)
        self.node_id = node_id
        self.gene_number = node_id[1]
        # Convert from stupid attribute type
        self.data = dict([(k, v) for k, v in data.items()])
        self.gene = self.diagram.net.genes[self.gene_number]
        self.text = "%s-%s" % node_id


class TextGeneShape(Shape):
    width = 2.2
    height = .55

    def __init__(self, diag, px, py, node_id, data):
        Shape.__init__(self, diag, px, py)
        self.outline = rounded_rect(px, py, self.width, self.height)
        self.gene_number = node_id[1]
        self.gene = self.diagram.net.genes[self.gene_number]
        self.data = data
        # self.text = self.gene.description
        # self.text = "%d" % self.gene_number
        p = self.diagram.net.collection.parameters
        varnames = [signal_to_variable(p, s) for s in self.gene.sub]
        varnames.reverse()  # This is just the way things are
        txt = self.gene.raw_description.format(*varnames)
        txt = "%s|%s" % (txt, signal_to_variable(p, self.gene.pub))
        self.text = latexify(txt)
        # latexify(self.gene.description)


class BaseDiagram(Diagram):
    def __init__(self, color=True):
        # self.scheme = ColorScheme if color else GrayScheme
        self.scheme = GrayScheme

    def draw(self, C):
        self.layout.draw(C)


class AChannelShape(Shape):
    size = .55
    radius = .35

    def __init__(self, diag, px, py, node_id):
        Shape.__init__(self, diag, px, py)
        self.node_id = node_id
        self.channel = node_id[1]
        if diag.graph.is_input(self.node_id):
            self.outline = cut_hexagon(px, py, self.radius)
            self.color = ColorScheme.cue
        elif diag.graph.is_output(self.node_id):
            self.outline = ppath.circle(px, py, self.radius)
            self.color = ColorScheme.act
        else:
            self.outline = diamond(px, py, self.radius)
            self.color = ColorScheme.sig

        self.text = "%s-%s" % node_id

    def draw(self, c):
        c.fill(self.outline, [self.color])
        self.draw_text(c, reversetext=True)


class ChannelShape(Shape):
    size = .55
    radius = .35
    diamond_radius = .55

    def __init__(self, diag, px, py, node_id):
        Shape.__init__(self, diag, px, py)
        # self.outline = diamond(px, py, self.size)
        # self.outline = hexagon(px, py, self.radius)
        self.node_id = node_id
        self.channel = node_id[1]
        if True:
            self.outline = diamond(px, py, self.diamond_radius)
        else:
            if diag.graph.is_input(self.node_id):
                self.outline = cut_hexagon(px, py, self.radius)
                self.color = ColorScheme.cue
            elif diag.graph.is_output(self.node_id):
                self.outline = ppath.circle(px, py, self.radius)
                self.color = ColorScheme.act
            else:
                self.outline = hexagon(px, py, self.radius)
                self.color = ColorScheme.sig

        # p = self.diagram.net.collection.parameters
        # c = self.channel
        self.text = 'test'
        # self.text = signal_to_variable(p, c)
        # self.text = "%s" % p.char_for_signal[self.channel]
        # self.text = "%s$_%s$" % (p.char_for_signal[self.channel], signal)

    def draw(self, c):
        if True:
            c.stroke(self.outline)  # , [pcolor.cmyk.White])
            self.draw_text(c)
        else:
            c.fill(self.outline, [self.color])
            self.draw_text(c, reversetext=True)

            # on = False
            # if self.diagram._phase == PhaseType.BIND:
            # on = elems[self.channel]
            # if not on:
            #     # Check out what's happening in the genes
            #     for g in self.incoming:
            #         if self.diagram._active[g.gene_number]:
            #             on = True
            #             break

            # Diamond.draw(self, c, reversetext=on)


class GeneConnector(CurvePath):
    # Map the binding types to arrows
    #
    def __init__(self, s1, s2, points, data):
        CurvePath.__init__(self, s1, s2, points)
        # TODO: Just defaulting to arrows for now
        self.edge_data = data
        # self.deco_fun = self.arrow_types[data['kind']]
        self.deco_fun = CurvePath.arrow_end

        if 'new' in data:
            col = pcolor.cmyk.Red
        elif 'old' in data:
            col = pcolor.cmyk.Gray
        else:
            col = pcolor.rgb.black
        self.color = col

    def draw(self, C):
        active = False
        # if self.diagram._phase == PhaseType.BIND:
        #     if isinstance(self.shape1, ChannelShape):
        #         if self.diagram._elements[self.shape1.channel]:
        #             active = True
        # else:
        #     if isinstance(self.shape2, ChannelShape):
        #         if self.diagram._active[self.shape1.gene_number]:
        #             active = True

        CurvePath.draw(self, C, active, self.color)


class BaseDotLayout(object):
    def write(self, f):
        self.agraph.write(f)

    def draw(self, f):
        self.agraph.draw(f, prog='dot', args=_dot_default_args)

    def get_gene_description(self, g):
        # return "G%02d: %s" % (g.index, g.description)
        return "%s" % (g.description)


class SimpleLayout(BaseLayout):
    def __init__(self, diagram):
        BaseLayout.__init__(self, diagram)
        # Compress along the yscale a bit
        self.yscaling = .015
        self.xscaling = .010

    def get_gene_description(self, g):
        # This is key for how the layout looks ... it gives space
        return "X"

    def get_gene_shape(self, node_id, px, py, data):
        return SimpleGeneShape(self.diagram, px, py, node_id, data)

    def get_signal_shape(self, node_id, px, py, data):
        return ChannelShape(self.diagram, px, py, node_id)

    def get_connection(self, shape1, shape2, points, data):
        return GeneConnector(shape1, shape2, points, data)


class TextLayout(BaseLayout):
    def __init__(self, diagram):
        BaseLayout.__init__(self, diagram)
        # Compress along the yscale a bit

    def get_gene_shape(self, node_id, px, py, data):
        return TextGeneShape(self.diagram, px, py, node_id, data)

    def get_signal_shape(self, node_id, px, py, data):
        return ChannelShape(self.diagram, px, py, node_id)

    def get_connection(self, shape1, shape2, points, data):
        return GeneConnector(shape1, shape2, points, data)


class AnalysisLayout(BaseLayout):
    def __init__(self, diagram):
        BaseLayout.__init__(self, diagram)
        # Compress along the yscale a bit
        self.yscaling = .015
        self.xscaling = .025

    def get_gene_description(self, g):
        # This is key for how the layout looks ... it gives space
        return "X-1"

    def get_gene_shape(self, node_id, px, py, data):
        return AGeneShape(self.diagram, px, py, node_id, data)

    def get_signal_shape(self, node_id, px, py, data):
        return AChannelShape(self.diagram, px, py, node_id)

    def get_connection(self, shape1, shape2, points, data):
        return GeneConnector(shape1, shape2, points, data)


class WiringDiagram(BaseDiagram):
    def __init__(self, net, color=True, layout_type=SimpleLayout):
        BaseDiagram.__init__(self, color)
        self.net = net
        self.graph = WiringGraph(net)
        self.graph.clean_graph()

        self.layout = layout_type(self)
        dm = DotMaker(self.graph)
        dm.make_layout(self.layout)


class SignallingDiagram(BaseDiagram):
    def __init__(self, net, color=True, layout_type=SimpleLayout):
        BaseDiagram.__init__(self, color)
        self.net = net
        self.graph = SignalGraph(net)
        self.graph.clean_graph()

        self.layout = layout_type(self)
        dm = DotMaker(self.graph)
        dm.make_layout(self.layout)


def test():
    import logging
    logging.basicConfig()
    nets = networks_from_yaml_file('tests/networks.yaml')
    # for i, n in enumerate(nets):
    #     d = SignallingDiagram(n)
    #     C = canvas.canvas()
    #     d.draw_basic(C)
    #     C.writePDFfile("net-%02d" % i)
    #
    # d = SignallingDiffDiagram(n1, n2)
    n1, n2 = nets[7], nets[8]
    n2 = nets[3]
    C = canvas.canvas()
    d = SignallingDiagram(n2, layout_type=TextLayout)
    d.draw(C)
    C.writePDFfile('signal')

    C = canvas.canvas()
    d = SignallingDiagram(n2, layout_type=SimpleLayout)
    d.draw(C)
    C.writePDFfile('signal2')

    d = WiringDiagram(n2, layout_type=TextLayout)
    C = canvas.canvas()
    d.draw(C)
    C.writePDFfile('wire')
    # st = d.get_states()
    # for s in st:
    #     s.draw(C)
    # C.writePDFfile('step')
    #


if __name__ == "__main__":
    test()
