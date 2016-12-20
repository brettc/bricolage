import numpy
from bricolage.analysis import NetworkAnalysis, AverageControlAnalyzer
from bricolage.graph_maker import encode_module_id, decode_module_id, NodeType
# from bricolage import dot_layout


def get_max_bindings(net):
    p = 0
    for g in net.genes:
        for m in g.modules:
            p += m.site_count()
    return p


def test_bindings(bowtie_network):
    net = bowtie_network
    a = NetworkAnalysis(net)
    a.get_edges()
    a.get_active_edges()
    max_bindings = get_max_bindings(net)
    assert a.active_bindings != 0
    assert a.potential_bindings <= max_bindings
    assert a.active_bindings <= a.potential_bindings


def test_bindings_pop(bowtie_database, bowtie_network):
    pop = bowtie_database.population

    max_bindings = get_max_bindings(bowtie_network)

    all_bindings = pop.active_bindings
    for i, net in enumerate(pop):
        ana = NetworkAnalysis(net)
        ana.get_active_edges()
        ana.get_edges()
        assert all_bindings[i] == ana.active_bindings

        # assert ana.active_bindings != 0
        # One of these turns out to be zero. Wow!
        # if ana.active_bindings == 0:
        #     print ana.get_active_edges()
        #     print ana.get_edges()
        #     dot_layout.save_network_as_fullgraph(net, name='no-bind', simplify=False)
        #     print net.attractors
        #     print net.rates

        assert ana.potential_bindings <= max_bindings
        assert ana.active_bindings <= ana.potential_bindings


def test_analysis_control(bowtie_database):
    az = AverageControlAnalyzer(bowtie_database.world)
    info = az.calc_info(bowtie_database.population)

    # Try accessing these
    info.control.mean(axis=0)
    info.entropy.mean(axis=0)
    print info.with_control


def test_module_id_encoding():
    for i in range(10):
        for j in range(3):
            x = encode_module_id(i, j)
            di, dj = decode_module_id(x)
            assert i == di
            assert j == dj


def test_active_cis(bowtie_database):
    pop = bowtie_database.population
    c_arr = pop.active_cis()

    cis_size = pop.factory.gene_count * pop.factory.module_count
    p_arr = numpy.zeros(cis_size, float)

    for net in pop:
        ana = NetworkAnalysis(net)
        for (et1, en1), (et2, en2) in ana.get_active_edges():
            if et2 == NodeType.MODULE:
                g, c = decode_module_id(en2)
                p_arr[g * pop.factory.module_count + c] += 1.0

    p_arr /= pop.size

    numpy.testing.assert_allclose(p_arr, c_arr)

