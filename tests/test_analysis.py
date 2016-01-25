from bricolage.analysis import NetworkAnalysis, AverageControlAnalyzer


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
        assert ana.active_bindings != 0
        assert ana.potential_bindings <= max_bindings
        assert ana.active_bindings <= ana.potential_bindings


def test_analysis_control(bowtie_database):
    az = AverageControlAnalyzer(bowtie_database.world)
    info = az.calc_info(bowtie_database.population)

    # Try accessing these
    info.control.mean(axis=0)
    info.entropy.mean(axis=0)
    print info.with_control

