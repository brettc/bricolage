import logging
logging.basicConfig()
import sys

from setpath import output_path, analysis_path
from invoke import task
from bricolage import threshold3, lineage, graph

def f1(a, b, c):
    return (a and not c) or (b and c)

def f2(a, b, c):
    return (a and not b) or (b and not c)

def make_func(a, b, c, d, e, f):
    r1 = [0, 0, 0]
    r2 = [0, 0, 0]
    if f1(a, b, c):
        r1 = [1, .5, .25]
    if f2(c, d, e):
        r2 = [1, .5, .25]
    return r1 + r2

def get_params(name, desc):
    return threshold3.Parameters(
        name=name,
        description=desc,
        seed=11, 
        cue_channels=6, 
        reg_channels=30,
        out_channels=6,
        cis_count=3,
        population_size=2000,
        mutation_rate=.001,
        replicates=20,
    )

def _create(rep, max_gen):
    lin = rep.get_lineage()
    if rep.fresh:
        lin.add_target(make_func)
    else:
        # Do some checking
        pass

    print '---> begin: ', rep.path
    print '** lineage is ', lin.path

    while lin.generation < max_gen:
        lin.next_generation()
        minf, maxf = lin.population.worst_and_best()
        if lin.generation % 25 == 0:
            sys.stdout.write("({}:{})".format(lin.generation, maxf))
            sys.stdout.flush()
        if maxf == 1.0:
            print
            raise StopReplicates
    print

    # When we're done, draw
    _draw(rep, lin, 5)
    print '---> end: ', rep.path

@task
def create(overwrite=False):
    name = "double"
    desc = "double"
    treat = lineage.Treatment(output_path / name, 
                              params=get_params(name, desc), 
                              analysis_path = analysis_path / name,
                              overwrite=overwrite,
                              full=False)
    treat.run(_create, max_gen=250000)

def _draw(rep, lin, maxn):
    winners=[]
    for net in lin.population:
        if net.fitness == 1.0:
            winners.append(net)

    for i, net in enumerate(winners):
        ana = logic2.NetworkAnalysis(net)
        g = graph.FullGraph(ana)#, knockouts=False)
        d = graph.DotMaker(g)
        p = str(rep.analysis_path / './net-{:02d}.png'.format(i))
        print '   ', p
        d.save_picture(p)
        if i > maxn:
            break

def _measure_mincut(rep, tfun, maxn=5):
    lin = rep.get_lineage()
    winners=[]
    for net in lin.population:
        if net.fitness == 1.0:
            winners.append(net)
        print net

@task 
def mincut_all():
    for i in range(len(get_functions())):
        name, tfun = get_functions()[i]
        name = "thresh-{:02d}".format(i)
        tfun = make_target_fun(tfun)
        try:
            treat = lineage.Treatment(output_path / name, params=get_params(name))
            treat.run(_measure_mincut, tfun)
                                    
        except:
            print "cannot find lineage {}".format(name)


def _calc_times(rep, gens, reps):
    lin = rep.get_lineage(readonly=True)
    g = lin._generations
    for i, row in enumerate(g.iterrows()):
        if row['best'] == 1.0:
            break
    else:
        i = None
    print lin.path, i
    gens.append(i)
    reps.append(rep.sequence)

