from math import sin, cos, pi
from pyx import trafo, deformer, canvas, path, deco, color, style, box
from pyx import text as pyx_text
from pyx.metapost.path import (beginknot, endknot, smoothknot, tensioncurve,
                               path as mpath)

import logtools

log = logtools.get_logger()

# We'll just set this stuff globally. Maybe it will matter later, but it
# doesn't right now.
# TODO: This appears to generate some annoying temp files
pyx_text.set(mode="latex")

# TODO: Make a pallette class here, that synchronises the colors with latex
# col = color.cmyk.PineGreen
# text.preamble(r"\usepackage{color}")
# text.preamble(r"\definecolor{COL}{cmyk}{%g,%g,%g,%g}" % (col.c, col.m,
# col.y, col.k))
# c = canvas.canvas()
# c.text(0, 0, r"\textcolor{COL}{Text} and outline have the same color")
# cmyk.Black          = cmyk(0, 0, 0, 1)
# cmyk.White          = cmyk(0, 0, 0, 0)


def make_preamble(runner):
    # runner.preamble(r"\usepackage{cmbright}")
    # We use the Helvetica font as it is everywhere (we can import into
    # Illustrator etc)
    runner.preamble(
        r"\usepackage[T1]{fontenc}\usepackage[scaled]{"
        r"helvet}\renewcommand*\familydefault{\sfdefault}\usepackage{sfmath}")
    runner.preamble(r"\parindent=0pt")
    runner.preamble(r"\usepackage{color}")
    runner.preamble(r"""
    \definecolor{white}{rgb}{1, 1, 1}
    \definecolor{black}{rgb}{0, 0, 0}
    \definecolor{blue}{RGB}{60, 147, 248}
    \definecolor{red}{RGB}{227, 69, 69}
    """)


make_preamble(pyx_text)

# TODO: Make this into a class, add color stuff to it
_runner = None
_cache = {}


def get_text(text):
    global _runner, _cache
    if _runner is None:
        _runner = pyx_text.texrunner(mode='latex')
        make_preamble(_runner)

    try:
        result = _cache[text]
    except KeyError:
        result = _runner.text(0, 0, text)
        _cache[text] = result

    return result


def rounded_rect(px, py, w, h, radius=.1):
    p = path.rect(px - w / 2.0, py - h / 2.0, w, h)
    s = deformer.cornersmoothed(radius)
    p = s.deform(p)
    return p


def diamond(px, py, w):
    p = path.rect(px - w / 2.0, py - w / 2.0, w, w)
    return p.transformed(trafo.rotate(45, px, py))


def hexagon(px, py, r):
    n = 6
    p = box.polygon(
        [(-r * sin(i * 2 * pi / n), r * cos(i * 2 * pi / n))
         for i in range(n)])
    p.transform(trafo.translate(px, py))
    return p.path()


def cut_hexagon(px, py, r):
    def lift(n):
        x, y = pts[n]
        y += r / 5.0
        pts[n] = x, y

    n = 6
    pts = [(-r * sin(i * 2 * pi / n), r * cos(i * 2 * pi / n))
           for i in range(n)]
    del pts[0:1]
    lift(0)
    lift(-1)
    p = box.polygon(pts)
    p.transform(trafo.translate(px, py))
    return p.path()


class Shape(object):
    alignment = [pyx_text.valign.middle, pyx_text.halign.center]

    def __init__(self, diag, px, py, text=None):
        self.diagram = diag
        self.px = px
        self.py = py
        self.text = text
        self.outgoing = []
        self.incoming = []

    def draw(self, c, reversetext=False):
        if reversetext:
            c.stroke(self.outline, [deco.filled([color.rgb.black])])
        else:
            c.stroke(self.outline, [deco.filled([color.rgb.white])])

        if self.text:
            self.draw_text(c, reversetext)

    def draw_text(self, c, reversetext=False):
        if reversetext:
            t = r"\textcolor{white}{%s}" % self.text
        else:
            t = self.text
        c.text(self.px, self.py, t, self.alignment)


class RoundedRect(Shape):
    def __init__(self, diag, px, py, w, h, text=None):
        Shape.__init__(self, diag, px, py, text)
        self.outline = rounded_rect(px, py, w, h)


class Diamond(Shape):
    def __init__(self, diag, px, py, w, text=None):
        Shape.__init__(self, diag, px, py, text)
        self.outline = diamond(px, py, w)


class Circle(Shape):
    def __init__(self, diag, px, py, r, text=None):
        Shape.__init__(self, diag, px, py, text)
        self.outline = path.circle(px, py, r)


class Connector(object):
    t_width = .12
    arrow_arc = t_width * 1.8
    arrow_angle = 18.0

    def __init__(self, s1, s2, deco_fun):
        self.diagram = s1.diagram
        self.shape1 = s1
        self.shape2 = s2
        self.deco_fun = deco_fun

        s1.outgoing.append(s2)
        s2.incoming.append(s1)

    def tee_end(self, pyx_canvas, options):
        # arc length and coordinates of tip
        ex, ey = self.arc.atend()
        bx, by = self.arc.at(self.arc.arclen() - self.t_width)
        line = path.line(bx, by, ex, ey)

        # from this template, we construct the two outer curves of the arrow
        x1 = line.transformed(trafo.rotate(-90, ex, ey))
        x2 = line.transformed(trafo.rotate(90, ex, ey))
        cross = x1.reversed() << x2

        pyx_canvas.stroke(cross, options + [style.linewidth.Thick])

    def arrow_end(self, pyx_canvas, options):
        # arc length and coordinates of tip
        ex, ey = self.arc.atend()
        bx, by = self.arc.at(self.arc.arclen() - self.arrow_arc)
        line = path.line(bx, by, ex, ey)

        # from this template, we construct the two outer curves of the arrow
        x1 = line.transformed(trafo.rotate(-self.arrow_angle, ex, ey))
        x2 = line.transformed(trafo.rotate(self.arrow_angle, ex, ey))
        arrow = x1.reversed() << x2
        arrow[-1].close()
        pyx_canvas.fill(arrow, options)
        pyx_canvas.stroke(arrow, options)
        # + [deco.filled([color.rgb.black])])

    def clip_to_nodes(self):
        p = self.arc

        if self.shape1 is self.shape2:
            # Must be a loop
            _, f_intersect = self.shape1.outline.intersect(p)
            _, p, _ = p.split(f_intersect)
            p, _ = p.split([p.end() - 0.1])
        else:

            # Clip the beginning, and clip the end
            _, f_intersect = self.shape1.outline.intersect(p)
            _, p = p.split(f_intersect)

            _, t_intersect = self.shape2.outline.intersect(p)
            p, _ = p.split(t_intersect)

            # Now clip it even more so that there is a gap on join
            p, _ = p.split([p.end() - 0.1])

        self.arc = p

    def draw(self, pyx_canvas, active=False, color=color.rgb.black):
        options = [color]
        if active:
            options.append(style.linewidth.THick)

        pyx_canvas.stroke(self.arc, options)
        if self.deco_fun:
            self.deco_fun(self, pyx_canvas, options)


class LinePath(Connector):
    def __init__(self, s1, s2, deco_fun=None):
        Connector.__init__(self, s1, s2, deco_fun)
        self.arc = path.line(s1.px, s1.py, s2.px, s2.py)
        self.clip_to_nodes()


class CurvePath(Connector):
    def __init__(self, s1, s2, points, deco_fun=None):
        Connector.__init__(self, s1, s2, deco_fun)

        knots = [beginknot(s1.px, s1.py)]
        knots.append(tensioncurve(ltension=2.0))
        for p in points:
            knots.append(smoothknot(*p))
            knots.append(tensioncurve())
        knots.append(endknot(s2.px, s2.py))

        # Make a metapost path from this
        self.arc = mpath(knots)
        self.clip_to_nodes()


class Diagram(object):
    def __init__(self):
        # self.smoother = deformer.cornersmoothed(.1)
        self._draw_list = []

    def add_object(self, obj, zorder=1):
        self._draw_list.append((zorder, obj))

    def draw(self, pyx_canvas):
        self._draw_list.sort()
        for z, obj in self._draw_list:
            obj.draw(pyx_canvas)

    # def connect_line(self, f, t):
    #     return LinePath(f, t)
    #
    # def connect_curve(self, f, t, points):
    #     return CurvePath(f, t, points)


def test1():
    d = Diagram()
    a = Diamond(d, 1, 2, .8, 'A')
    b = Circle(d, 5, 6, .4, 'B')
    c = RoundedRect(d, 1, 6, 2, 1, r'$A_i$ $\land{}$ A$_j$')
    # d = HalfRoundedRect(d, 10, 8, 2, 1)

    pyx_canvas = canvas.canvas()
    a.draw(pyx_canvas)
    b.draw(pyx_canvas, True)
    c.draw(pyx_canvas, True)
    # d.draw(pyx_canvas)

    LinePath(a, b, Connector.arrow_end).draw(pyx_canvas)
    # LinePath(a, c).draw(pyx_canvas)
    # CurvePath(c, b, [(3, 3)], Connector.tee_end).draw(pyx_canvas)

    pyx_canvas.writePDFfile("testing")


if __name__ == '__main__':
    test1()
