import logging

log = logging.getLogger("drawing")

from graph_maker import ChannelType
from graph_layout import DotDiagram
from pyx_drawing import (Diagram, rounded_rect, CurvePath, Shape, cut_hexagon,
                       hexagon, diamond)
from pyx import canvas, color as pcolor, path as ppath
from .cis_logic import text_for_gene, text_for_cis_mod

def latexify(l):
    l = l.replace("~", r"$\lnot{}$")
    # l = l.replace("~", r"$\sim$")
    # l = l.replace("~", "NOT ")
    # l = l.replace("~", r"\~{}")
    l = l.replace(" & ", r" $\land{}$ ")
    l = l.replace(" | ", r" $\lor{}$ ")
    # l = l.replace(")", "")
    # l = l.replace("(", "")
    l = l.replace("|", r"$\to$")
    return l


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


class SimpleGeneShape(Shape):
    width = .8
    height = .2

    def __init__(self, diag, px, py, gene):
        Shape.__init__(self, diag, px, py)
        self.outline = rounded_rect(px, py, self.width, self.height)

    def draw(self, pyx_canvas):
        pyx_canvas.fill(self.outline, [pcolor.cmyk.Gray])


class SimpleChannelShape(Shape):
    size = .55
    radius = .35

    def __init__(self, diag, px, py, channel, channel_type):
        Shape.__init__(self, diag, px, py)
        self.channel = channel
        if channel_type == ChannelType.INPUT:
            self.outline = cut_hexagon(px, py, self.radius * .8)
            self.color = ColorScheme.cue
        elif channel_type == ChannelType.OUTPUT:
            self.outline = ppath.circle(px, py, self.radius * .8)
            self.color = ColorScheme.act
        else:
            self.outline = diamond(px, py, self.radius)
            self.color = ColorScheme.sig
        self.text = str(channel)

    def draw(self, c):
        c.fill(self.outline, [self.color])
        # self.draw_text(c, reversetext=True)


class TextGeneShape(Shape):
    width = 2.8
    height = .55

    def __init__(self, diag, px, py, gene):
        Shape.__init__(self, diag, px, py)
        self.outline = rounded_rect(px, py, self.width, self.height)
        self.text = latexify(text_for_gene(gene.network.factory.world, gene))
        # latexify(self.gene.description)


class SimpleDiagram(DotDiagram):
    def __init__(self, graph):
        DotDiagram.__init__(self, graph, yscale=.8, xscale=.7)

    # TODO!!!
    # def get_gene_description(self, g):
    #     # This is key for how the layout looks ... it gives space
    #     return "X"

    def get_gene_shape(self, px, py, gene):
        return SimpleGeneShape(self, px, py, gene)

    def get_signal_shape(self, px, py, channel, channel_type):
        return SimpleChannelShape(self, px, py, channel, channel_type)

    def get_connector(self, shape1, shape2, points):
        return GeneConnector(shape1, shape2, points)




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
    def __init__(self, s1, s2, points):
        CurvePath.__init__(self, s1, s2, points)
        self.deco_fun = CurvePath.arrow_end
        self.color = pcolor.rgb.black

    def draw(self, pyx_canvas):
        active = False
        CurvePath.draw(self, pyx_canvas, active, self.color)


# class TextLayout(BaseLayout):
#     def __init__(self, diagram):
#         BaseLayout.__init__(self, diagram)
#         # Compress along the yscale a bit
#
#     def get_gene_shape(self, node_id, px, py, data):
#         return TextGeneShape(self.diagram, px, py, node_id, data)
#
#     def get_signal_shape(self, node_id, px, py, data):
#         return ChannelShape(self.diagram, px, py, node_id)
#
#     def get_connector(self, shape1, shape2, points, data):
#         return GeneConnector(shape1, shape2, points, data)
#
#
# class AnalysisLayout(BaseLayout):
#     def __init__(self, diagram):
#         BaseLayout.__init__(self, diagram)
#         # Compress along the yscale a bit
#         self.yscaling = .015
#         self.xscaling = .025
#
#     def get_gene_description(self, g):
#         # This is key for how the layout looks ... it gives space
#         return "X-1"
#
#     def get_gene_shape(self, node_id, px, py, data):
#         return AGeneShape(self.diagram, px, py, node_id, data)
#
#     def get_signal_shape(self, node_id, px, py, data):
#         return AChannelShape(self.diagram, px, py, node_id)
#
#     def get_connector(self, shape1, shape2, points, data):
#         return GeneConnector(shape1, shape2, points, data)
#

def test():
    import logging
    logging.basicConfig()
    nets = networks_from_yaml_file('tests/networks.yaml')
    # for i, n in enumerate(nets):
    #     d = SignallingDiagram(n)
    #     pyx_canvas = canvas.canvas()
    #     d.draw_basic(pyx_canvas)
    #     pyx_canvas.writePDFfile("net-%02d" % i)
    #
    # d = SignallingDiffDiagram(n1, n2)
    n1, n2 = nets[7], nets[8]
    n2 = nets[3]
    pyx_canvas = canvas.canvas()
    d = SignallingDiagram(n2, layout_type=TextLayout)
    d.draw(pyx_canvas)
    pyx_canvas.writePDFfile('signal')

    pyx_canvas = canvas.canvas()
    d = SignallingDiagram(n2, layout_type=SimpleLayout)
    d.draw(pyx_canvas)
    pyx_canvas.writePDFfile('signal2')

    d = WiringDiagram(n2, layout_type=TextLayout)
    pyx_canvas = canvas.canvas()
    d.draw(pyx_canvas)
    pyx_canvas.writePDFfile('wire')
    # st = d.get_states()
    # for s in st:
    #     s.draw(pyx_canvas)
    # pyx_canvas.writePDFfile('step')
    #


if __name__ == "__main__":
    test()
