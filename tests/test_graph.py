from bricolage import graph, threshold3
from bricolage.logic_tools import text_for_gene, text_for_cis_mod
import pathlib
from bricolage.core_ext import Network
from bricolage.analysis_ext import NetworkAnalysis
from bricolage.core import InterventionState


def test_graph_creation(tmpdir, bowtie_network):
    output_path = pathlib.Path(str(tmpdir))
    ana = NetworkAnalysis(bowtie_network)
    g = graph.SignalFlowGraph(ana)
    d = graph.DotMaker(g)
    d.save_picture(str(output_path / 'signal.png'))

    # g = graph.GeneSignalGraph(ana)
    g = graph.GeneGraph(ana)
    d = graph.DotMaker(g)
    d.save_picture(str(output_path / 'gene.png'))

    g = graph.FullGraph(ana)
    d = graph.DotMaker(g)
    d.save_picture(str(output_path / 'full.png'))


def test_1(bowtie_network):
    net = bowtie_network
    w = net.factory.world
    assert isinstance(net, Network)
    ana = NetworkAnalysis(net)
    ana.get_active_edges()
    mod = ana.modified
    print

    g = graph.GeneSignalGraph(ana, knockouts=False)
    d = graph.DotMaker(g)
    d.save_picture('./full.png')

    g = graph.GeneSignalGraph(ana)
    d = graph.DotMaker(g)
    d.save_picture('./full_x.png')

    for go, gm in zip(net.genes, mod.genes):
        if InterventionState.INTERVENE_OFF == gm.intervene:
            continue

        print text_for_gene(w, go), '---', text_for_gene(w, gm)

        # for mo, mm in zip(go.modules, gm.modules):
        #     if InterventionState.INTERVENE_OFF == mm.intervene:
        #         continue
        #     # print mo.intervene, mm.intervene
        #     # for co, cm in zip(mo.channels, gm.channels):
        #     print go,
        #     # print mo.channels, mm.channels
        #     compute_boolean(w, mm)


def test_boolean_binding():
    params = threshold3.Parameters(
        cis_count=3,
        reg_channels=6,
        out_channels=3,
        cue_channels=3,
        population_size=100,
        mutation_rate=.001,
    )

    world = threshold3.World(params)
    print boolean_func_from_coop_binding(world, [2, 1, 4], [2, 2, 1])




