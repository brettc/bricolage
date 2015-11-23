from bricolage import analysis_ext, graph, threshold3
from bricolage.logic_tools import boolean_func_from_coop_binding
from generate import get_database
import pathlib

def test_graph_creation(tmpdir):
    db = get_database('bowtie', readonly=True)
    interesting_network = db.population.get_best()[0]

    output_path = pathlib.Path(str(tmpdir))
    ana = analysis_ext.NetworkAnalysis(interesting_network)
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
    db.close()

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

