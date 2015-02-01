import pytest
from bricolage import threshold3, lineage, graph
import pathlib
# import networkx

@pytest.fixture
def target_2x2():
    def f(a, b):
        if a and b:
            x = .5
        else:
            x = 0

        if a == b:
            y = 1.0
        else:
            y = 0

        return [x, y]
    return f

@pytest.fixture
def p_2x2():
    return threshold3.Parameters(
        seed=5, 
        cis_count=3, 
        reg_channels=5,
        out_channels=2, 
        cue_channels=2, 
        population_size=1000,
        mutation_rate=.001,
    )

@pytest.fixture
def interesting_network(tmpdir, p_2x2, target_2x2):
    base = pathlib.Path(str(tmpdir))
    path = base / 'signal_flow.db'
    a = lineage.SnapshotLineage(path, params=p_2x2)
    a.add_target(target_2x2)

    # Get some nice graphs
    for i in range(100):
        a.next_generation()
        w, b = a.population.worst_and_best()
        if b == 1.0:
            break

    for n in a.population:
        if n.fitness == 1.0:
            break
    
    return n

def test_graph_creation(tmpdir, interesting_network):
    output_path = pathlib.Path(str(tmpdir))
    ana = threshold3.NetworkAnalysis(interesting_network)
    g = graph.SignalFlowGraph(ana)
    for n in g.minimum_cut():
        print g.get_channel_label(n[1])
    d = graph.DotMaker(g)
    d.save_picture(str(output_path / 'signal.png'))

    # g = graph.GeneSignalGraph(ana)
    g = graph.GeneGraph(ana)
    d = graph.DotMaker(g)
    d.save_picture(str(output_path / 'gene.png'))

    g = graph.FullGraph(ana)
    d = graph.DotMaker(g)
    d.save_picture(str(output_path / 'full.png'))


