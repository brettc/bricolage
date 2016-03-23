import numpy
from generate import make_noisy, get_numpy_dump

def test_noisy_fitnesses():
    f1 = make_noisy()
    f2 = get_numpy_dump('noisy')
    numpy.testing.assert_allclose(f1, f2)


