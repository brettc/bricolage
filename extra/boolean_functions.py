import itertools
from sympy.logic import SOPform, bool_map
from sympy.printing.latex import latex
# sympy.logic.boolalg.bool_map(bool1, bool2)[source]

def detex(term):
    rt = latex(term)
    rt = rt.replace(r"\left", "")
    rt = rt.replace(r"\right", "")
    rt = rt.replace(r"\wedge", " and ")
    rt = rt.replace(r"\vee", " or ")
    rt = rt.replace(r"\neg ", " not ")
    rt = " ".join(rt.strip().split())
    # rt = rt.replace(" ", "")
    return rt

def min_terms_for(inputs, op):
    mt = []
    for inp, out in zip(inputs, op):
        if out:
            mt.append(inp)
    return mt

# This prunes out any structures which are "equivalent" (ie. swapping a/b/c
# around makes exactly the same statement
def already_in_here(soplist, sop):
    for sopother in soplist:
        sameas = bool_map(sop, sopother)
        # Not sure why we get True / False or tuple. Weird?
        if isinstance(sameas, tuple):
            print 'these are the same', detex(sop), '---', detex(sopother)
            return True
    return False

# Produce all of the basic 3 place boolean functions where half the function
# is True, the rest False.
#
# Prune out the ones which don't rely on all three inputs
def function_maker():
    sop_list = []
    func_list = []
    lookup_list = []

    # Every possible output function for 3 input bool func
    all_poss = itertools.product([0, 1], repeat=8)

    # All possible inputs
    inputs = [_ for _ in itertools.product([0, 1], repeat=3)]

    for resp in all_poss:
        # Only look at those with equal True/False (as this gives us a nice
        # information measure)
        if sum(resp) == 4:
            # Make sure it uses all the variables
            mt = min_terms_for(inputs, resp)
            print mt
            sop = SOPform('a b c'.split(), mt)
            if len(sop.atoms()) == 3:
                # We found one that uses all 3
                if not already_in_here(sop_list, sop):
                    lookup = dict(zip(inputs, resp))
                    # def the_func(*ip):
                    #     c = lookup.copy()
                    #     return c[ip]
                    # fname = "f" + "".join([str(_) for _ in resp])
                    # the_func.__doc__ = detex(sop)
                    func_list.append(detex(sop))
                    sop_list.append(sop)
                    lookup_list.append(lookup)

    return func_list, lookup_list


if __name__ == '__main__':
    import cPickle as pickle
    names, function_map = function_maker()
    inputs = [_ for _ in itertools.product([0, 1], repeat=3)]
    for i, (f, l) in enumerate(zip(names, function_map)):
        print i, f
        # for j in inputs:
        #     print j, l[j]
        print '--'

    print 'saving to file'
    with open('minfuncs.pickle', 'wb') as minfuncs:
        pickle.dump([names, function_map], minfuncs, -1)



#
