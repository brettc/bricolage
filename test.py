from tests.generate import get_database
from bricolage.graph_layout import DotMaker
from bricolage.graph_draw import SimpleDiagram
from bricolage import graph_maker
from bricolage.analysis import NetworkAnalysis
from pyx import canvas

with get_database('bowtie', readonly=True) as db:
    net = db.population.get_best()[0]
    ana = NetworkAnalysis(net)

    t = db.targets[0]
    ana.calc_mutual_info(t)
    ana.calc_output_control()
    grph = graph_maker.GeneSignalGraph(ana)
    DotMaker(grph).save_picture('test.png')

    # grph = graph_maker.GeneSignalGraph(ana)
    # diag = SimpleDiagram(grph)
    # #
    # c = canvas.canvas()
    # diag.draw(c)
    # c.writePDFfile('test.pdf')
