import pathlib

from bricolage.dot_layout import DotMaker, DotDiagram
# from bricolage.graph_draw import SimpleLayout
from bricolage.pyx_drawing import Diagram
from bricolage import graph_maker, threshold3
from bricolage.cis_logic import text_for_gene
from bricolage.core_ext import Network
from bricolage.analysis_ext import NetworkAnalysis
from bricolage.core import InterventionState


def test_graph_creation(tmpdir, bowtie_network):
    output_path = pathlib.Path(str(tmpdir))
    ana = NetworkAnalysis(bowtie_network)
    g = graph_maker.SignalFlowGraph(ana)
    d = DotMaker(g)
    d.save_picture(str(output_path / 'signal.png'))

    # g = graph.GeneSignalGraph(ana)
    g = graph_maker.GeneGraph(ana)
    d = DotMaker(g)
    d.save_picture(str(output_path / 'gene.png'))

    g = graph_maker.FullGraph(ana)
    d = DotMaker(g)
    d.save_picture(str(output_path / 'full.png'))


# TODO: reinstate some tests
def test_layouts(bowtie_network):
    net = bowtie_network
    ana = NetworkAnalysis(net)
    grph = graph_maker.GeneSignalGraph(ana)
    dotdia = DotMaker(grph)
    dotm.make_diagram(SimpleLayout(Diagram()))


# def test_1(bowtie_network):
#     net = bowtie_network
#     w = net.factory.world
#     assert isinstance(net, Network)
#     ana = NetworkAnalysis(net)
#     ana.get_active_edges()
#     mod = ana.modified
#     print
#
#     g = graph_maker.GeneSignalGraph(ana, knockouts=False)
#     d = bricolage.graph_layout.DotMaker(g)
#     d.save_picture('./full.png')
#
#     g = graph_maker.GeneSignalGraph(ana)
#     d = bricolage.graph_layout.DotMaker(g)
#     d.save_picture('./full_x.png')
#
#     for go, gm in zip(net.genes, mod.genes):
#         if InterventionState.INTERVENE_OFF == gm.intervene:
#             continue
#
#         print text_for_gene(w, go), '---', text_for_gene(w, gm)
#
#         # for mo, mm in zip(go.modules, gm.modules):
#         #     if InterventionState.INTERVENE_OFF == mm.intervene:
#         #         continue
#         #     # print mo.intervene, mm.intervene
#         #     # for co, cm in zip(mo.channels, gm.channels):
#         #     print go,
#         #     # print mo.channels, mm.channels
#         #     compute_boolean(w, mm)

