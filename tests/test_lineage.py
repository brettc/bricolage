import bricolage.logic2 as T
import bricolage.lineage as L
import pathlib

def test_creation(tmpdir):
    # fname = pathlib.Path(str(tmpdir)) / 'testing.db'
    fname = 'testing.db'
    params = T.Parameters(population_size=10)
    lineage = L.Lineage(fname, params, T.Constructor)
    del lineage
    # lineage.close()
    l2 = L.Lineage(fname)
    # print l2
    for net in l2.population:
        print net.identifier

    # lineage.close()

