import numpy
import pytest

from bricolage.analysis import NetworkAnalysis, AverageControlAnalyzer
from bricolage.graph_maker import encode_module_id, decode_module_id, NodeType


def get_max_bindings(net):
    p = 0
    for g in net.genes:
        for m in g.modules:
            p += m.site_count()
    return p


def test_bindings(bowtie_network):
    net = bowtie_network
    a = NetworkAnalysis(net)
    # Force bindings to be calculated
    # TODO: this is total crap
    e = a.get_edges()
    ae = a.get_active_edges()
    max_bindings = get_max_bindings(net)
    assert a.active_bindings != 0
    assert a.potential_bindings <= max_bindings
    assert a.active_bindings <= a.potential_bindings


@pytest.mark.skip(reason="currently broken, not sure required")
def test_bindings_pop(bowtie_database, bowtie_network):
    pop = bowtie_database.population

    max_bindings = get_max_bindings(bowtie_network)
    all_bindings = pop.active_bindings

    # BRETT: Problematic here for some reason. The bindings are not the same.
    # Should the population one be the same. Or did I never get aronud to updateing the pop level ones.
    # wo
    # And things weirdly crash afterwards
    for i, net in enumerate(pop):
        print(i, net)
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


def test_module_id_encoding():
    for i in range(10):
        for j in range(3):
            x = encode_module_id(i, j)
            di, dj = decode_module_id(x)
            assert i == di
            assert j == dj


# Broken. Not used for now
@pytest.mark.skip(reason="currently broken, not sure required")
def test_active_cis(bowtie_database):
    pop = bowtie_database.population
    c_arr = pop.active_cis()
    print c_arr

    cis_size = pop.factory.gene_count * pop.factory.module_count
    p_arr = numpy.zeros(cis_size, float)

    for net in pop:
        ana = NetworkAnalysis(net)
        for (et1, en1), (et2, en2) in ana.get_active_edges():
            if et2 == NodeType.MODULE:
                g, c = decode_module_id(en2)
                p_arr[g * pop.factory.module_count + c] += 1.0
        break

    p_arr /= pop.size

    numpy.testing.assert_allclose(p_arr, c_arr)
