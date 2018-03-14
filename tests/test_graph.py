import pathlib
import pytest

from bricolage.dot_layout import DotMaker
from bricolage.graph_draw import SmallDiagram, TextDiagram
from bricolage import graph_maker
from bricolage.analysis_ext import NetworkAnalysis
from bricolage.core_ext import Collection
import numpy
import cPickle as pickle

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
    SmallDiagram(grph)
    TextDiagram(grph)

def test_save(bowtie_network, tmpdir):
    # tmpdir = pathlib.Path('.')

    with open('test.pickle', 'wb') as f:
        pickle.dump(bowtie_network, f, -1)

    with open('test.pickle', 'rb') as f:
        pickle.load(f)

    # print net2
    # print net2.factory.world
    #

def draw_net(pth, name, net):
    ana = NetworkAnalysis(net)
    g = graph_maker.GeneGraph(ana)
    d = DotMaker(g)
    d.save_picture((pth / name).with_suffix('.png').as_posix())


def graph_diff(net1, net2):
    ana1 = NetworkAnalysis(net1)
    g1 = graph_maker.FullGraph(ana1)
    ana2 = NetworkAnalysis(net2)
    g2 = graph_maker.FullGraph(ana2)

    # e1 = set([e for e in g1.nx_graph.edges()])
    # e2 = set([e for e in g2.nx_graph.edges()])
    # print e1.symmetric_difference(e2)
    #
    for node in g2.nx_graph.nodes():
        # Is it there?
        if node in g1.nx_graph.nodes():
            if g1.get_label(node).strip() != g2.get_label(node).strip():
                print 'changed', node
        else:
            print 'added', node

    for node in g1.nx_graph.nodes():
        if node not in g2.nx_graph.nodes():
            print 'removed', node

    for edge in g2.nx_graph.edges():
        if edge not in g1.nx_graph.edges():
            print 'added', edge

    for edge in g1.nx_graph.edges():
        if edge not in g2.nx_graph.edges():
            print 'removed', edge

    # n1 = set([g1.get_label(n).strip() for n in g1.nx_graph.nodes()])
    # n2 = set([g2.get_label(n).strip() for n in g2.nx_graph.nodes()])
    # print n1
    # print n1.symmetric_difference(n2)

@pytest.fixture
def net1(double_bow):
    return double_bow.population.get_best()[0]

def draw_diff(name, net1, net2):
    ana1 = NetworkAnalysis(net1)
    g1 = graph_maker.GeneGraph(ana1)
    ana2 = NetworkAnalysis(net2)
    g2 = graph_maker.GeneGraph(ana2)
    d = DotMaker(g2)
    a = d.get_dot_diff(g1)
    _dot_default_args = '-Nfontname="Helvetica-8"'
    a.draw(name, prog='dot', args=_dot_default_args)

def test_changes(net1):
    net = net1
    # output_path = pathlib.Path(str(tmpdir))
    output_path = pathlib.Path('.')
    c = Collection(net.factory)
    c.fill_with_mutations(net, numpy.asarray([1, 1, 2, 2]))
    other = c[1]

    draw_net(output_path, 'bowtie1', net)
    draw_net(output_path, 'bowtie2', other)
    graph_diff(net, other)

    draw_diff('diff1.png', net, other)
    draw_diff('diff2.png', other, net)

    # n1 = c[2]



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
#         print text_for_gene(go), '---', text_for_gene(gm)
#
#         # for mo, mm in zip(go.modules, gm.modules):
#         #     if InterventionState.INTERVENE_OFF == mm.intervene:
#         #         continue
#         #     # print mo.intervene, mm.intervene
#         #     # for co, cm in zip(mo.channels, gm.channels):
#         #     print go,
#         #     # print mo.channels, mm.channels
#         #     compute_boolean(w, mm)

