from bricolage import analysis_ext, graph
from generate import get_database
import pathlib

def test_graph_creation(tmpdir):
    db = get_database('bowtie')
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


