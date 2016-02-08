import logging

log = logging.getLogger("drawing")

from graph_maker import ChannelType, NodeType
from dot_layout import DotDiagram
from pyx_drawing import (rounded_rect, CurvePath, Shape, cut_hexagon,
                         cut_hexagon_rotated, measure_text)
from pyx import color as pcolor, path as ppath
from .cis_logic import text_for_gene
import re


RE_ENV = re.compile(r"E(\d+)")
RE_TRANS = re.compile(r"T(\d+)")
RE_STRUCT = re.compile(r"P(\d+)")


# EEEK. Truly awful.
def latexify(l):
    l = l.replace("~", r"\lnot{}")
    l = l.replace(" & ", r" \land{} ")
    l = l.replace(" | ", r" \lor{} ")
    l = l.replace("=>", r"\to{}")
    l = RE_ENV.sub(r"e_{\1}", l)
    l = RE_TRANS.sub(r"t_{\1}", l)
    l = RE_STRUCT.sub(r"p_{\1}", l)
    return "$" + l + "$"


class ColorScheme(object):
    sig = pcolor.cmyk.Red
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
        # txt = "e$_{%d}$" % (c - p.cue_range[0] + 1)
        txt = "%d" % (c - p.cue_range[0] + 1)
    elif c in p.reg_signals:
        # txt = "r$_{%d}$" % (c - p.reg_range[0] + 1)
        # txt = r"\textbf{%d}" % (c - p.reg_range[0] + 1)
        txt = r"%d" % (c - p.reg_range[0] + 1)
    else:
        # txt = "p$_{%d}$" % (c - p.out_range[0] + 1)
        txt = "%d" % (c - p.out_range[0] + 1)
    return txt


class GeneConnector(CurvePath):
    # Map the binding types to arrows
    def __init__(self, s1, s2, points):
        CurvePath.__init__(self, s1, s2, points)
        self.deco_fun = CurvePath.arrow_end
        self.color = pcolor.rgb.black

    def draw(self, pyx_canvas):
        active = False
        CurvePath.draw(self, pyx_canvas, active, self.color)


class SmallGeneShape(Shape):
    width = .8
    height = .2

    def __init__(self, diag, px, py, gene):
        Shape.__init__(self, diag, px, py)
        self.outline = rounded_rect(px, py, self.width, self.height)
        self.gene = gene
        if gene.sequence >= diag.world.reg_channels:
            self.color = pcolor.cmyk.Green
        else:
            self.color = pcolor.cmyk.Gray
            # self.outline = ppath.circle(px, py, .25)
            # self.text = str(gene.pub)

    def draw(self, pyx_canvas):
        pyx_canvas.fill(self.outline, [self.color])


class SmallChannelShape(Shape):
    size = .55
    radius = .35

    def __init__(self, diag, px, py, channel, channel_type):
        Shape.__init__(self, diag, px, py)
        self.channel = channel
        self.channel_type = channel_type

        if channel_type == ChannelType.INPUT:
            self.outline = cut_hexagon(px, py, self.radius * .8)
            self.color = ColorScheme.cue

        elif channel_type == ChannelType.OUTPUT:
            # self.outline = ppath.circle(px, py, self.radius * .8)
            self.outline = cut_hexagon_rotated(px, py, self.radius * .8)
            self.color = ColorScheme.act
            # self.color = pcolor.gradient.Jet.getcolor(channel * .05)
        else:
            self.outline = ppath.circle(px, py, self.radius * .8)
            # self.outline = diamond(px, py, self.radius)
            # self.color = pcolor.gradient.WhiteRed.getcolor(
            #     (channel - p.reg_range[0]) * (1.0/p.reg_channels))
            # self.color = pcolor.gradient.Gray.getcolor(annote[channel]['W'])
            self.color = pcolor.cmyk.Black
            # ColorScheme.sig

        if channel == diag.highlight:
            self.color = pcolor.cmyk.Red

        self.text = signal_to_variable(self.diagram.world, channel)

    def draw(self, c):
        c.fill(self.outline, [self.color])
        # if self.channel_type == ChannelType.INTERNAL:
        # c.stroke(self.outline, [pcolor.cmyk.Black])
        # self.draw_text(c)
        self.draw_text(c, reversetext=True)


class SmallDiagram(DotDiagram):
    def __init__(self, graph, height=9.0, width=5.0, highlight=None):
        self.highlight = highlight
        super(SmallDiagram, self).__init__(graph, height, width)

    def get_label(self, node_type, ident):
        # We just keep this as small as possible
        return 'X'

    def get_gene_shape(self, px, py, gene):
        return SmallGeneShape(self, px, py, gene)

    def get_signal_shape(self, px, py, channel, channel_type):
        return SmallChannelShape(self, px, py, channel, channel_type)

    def get_connector(self, shape1, shape2, points):
        return GeneConnector(shape1, shape2, points)


class TextGeneShape(Shape):
    height = .6

    def __init__(self, diag, px, py, gene):
        Shape.__init__(self, diag, px, py)
        self.text = latexify(
            text_for_gene(gene) + '=>' + diag.world.name_for_channel(gene.pub))
        w, h = measure_text(self.text)
        self.outline = rounded_rect(px, py, w + .2, self.height)
        if gene.sequence >= diag.world.reg_channels:
            self.color = pcolor.cmyk.Green
        else:
            self.color = pcolor.cmyk.Black

    def draw(self, c):
        c.stroke(self.outline, [self.color])
        self.draw_text(c)


class TextChannelShape(Shape):
    radius = .45

    def __init__(self, diag, px, py, channel, channel_type):
        Shape.__init__(self, diag, px, py)
        self.channel = channel
        self.channel_type = channel_type

        if channel_type == ChannelType.INPUT:
            self.outline = cut_hexagon(px, py, self.radius * .8)
            self.color = ColorScheme.cue

        elif channel_type == ChannelType.OUTPUT:
            # self.outline = ppath.circle(px, py, self.radius * .8)
            self.outline = cut_hexagon_rotated(px, py, self.radius * .8)
            self.color = ColorScheme.act
            # self.color = pcolor.gradient.Jet.getcolor(channel * .05)
        else:
            self.outline = ppath.circle(px, py, self.radius * .8)
            # self.outline = diamond(px, py, self.radius)
            # self.color = pcolor.gradient.WhiteRed.getcolor(
            #     (channel - p.reg_range[0]) * (1.0/p.reg_channels))
            # self.color = pcolor.gradient.Gray.getcolor(annote[channel]['W'])
            self.color = pcolor.cmyk.Black
            # ColorScheme.sig

        if channel == diag.highlight:
            self.color = pcolor.cmyk.Red

        self.text = latexify(self.diagram.world.name_for_channel(channel))

    def draw(self, c):
        c.fill(self.outline, [self.color])
        # if self.channel_type == ChannelType.INTERNAL:
        # c.stroke(self.outline, [pcolor.cmyk.Black])
        # self.draw_text(c)
        self.draw_text(c, reversetext=True)


class TextDiagram(DotDiagram):
    def __init__(self, graph, height=10.0, width=8.0, highlight=None):
        self.highlight = highlight
        super(TextDiagram, self).__init__(graph, height, width)

    def get_label(self, node_type, ident):
        # We just keep this as small as possible
        if node_type == NodeType.GENE:
            g = self.graph.network.genes[ident]
            return text_for_gene(g) + '=>' + self.world.name_for_channel(g.pub)

        return 'X'

    def get_gene_shape(self, px, py, gene):
        return TextGeneShape(self, px, py, gene)

    def get_signal_shape(self, px, py, channel, channel_type):
        return TextChannelShape(self, px, py, channel, channel_type)

    def get_connector(self, shape1, shape2, points):
        return GeneConnector(shape1, shape2, points)
