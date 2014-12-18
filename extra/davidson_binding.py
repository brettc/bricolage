"""
An explorating of how to convert cooperative binding strengths into boolean
operators. This generates (roughly) the kind of boolean operators that we see
in the Davidson-Peters models.
"""

import itertools
from sympy.logic import SOPform, simplify_logic
from sympy.printing.latex import latex

def term_to_simple_text(term):
    """This appears to be the easiest way to write the functions in a simple manner"""
    rt = latex(term)
    rt = rt.replace(r"\left", "")
    rt = rt.replace(r"\right", "")
    rt = rt.replace(r"\wedge", "&")
    rt = rt.replace(r"\vee", "|")
    rt = rt.replace(r"\neg ", "~")
    rt = rt.replace(r"\mathrm{", "")
    rt = rt.replace(r"}", "")
    rt = rt.replace(" ", "")
    rt = rt.replace(" ", "")
    return rt

def make_bindings(size, on_channel=False, thresh_adjust=0):
    """
    The rough idea. Each "site" either contributes to or detracts from the
    binding. Assuming there are three sites, We assign each site from -3 to +3.
    If a site is "bound", then we add the amount at that site. 
    """
    # Not sure this should be size or input size
    strengths = range(-size, size+1)
    # Assign each site a..z

    input_size = size
    if on_channel:
        input_size += 1
    letters = [chr(ord('a') + i) for i in range(input_size)]

    soplist = {}
    threshold = size - thresh_adjust

    # This generates all possible combinations of inputs for the number binary
    # sites
    inputs = list(itertools.product((0, 1), repeat=size))

    # This generates all possible combinations of binding strengths.
    for bindings in itertools.product(strengths, repeat=input_size):
        minterms = []

        # Iterate through all inputs for this
        for i in inputs:
            # If there is an on_channel (a background channel always on), then
            # we add this as channel 0
            if on_channel:
                cur_input = [1] + list(i)
            else:
                cur_input = i
            active = [a for a, b in zip(bindings, cur_input) if b == 1]
            if sum(active) >= threshold:
                minterms.append(i)

        # Now, what boolean function does that set of combinations work for?
        sop = SOPform(letters, minterms)

        # if len(sop.atoms()) > 0:
            # Some of them are repeated
        b = soplist.setdefault(sop, [])
        b.append(bindings)

    return soplist

def print_sops(sops, count=False):
    # Look at the distribution of boolean functions that are generated. Some occur
    # more than others ... 
    print "============="
    print "There are {} items!".format(len(sops))
    print "------------"
    if count:
        sort_sops = []
        for k, v in sops.iteritems():
            sort_sops.append((sum(v), k))
    else:
        sort_sops = [(len(v), k) for k, v in sops.items()]
    total = 0
    for number, eq in sorted(sort_sops, key=lambda x: x[0]):
        # print term_to_simple_text(eq), number
        print "${}$ | {}".format(latex(eq), number)
        total += number
    print "------------"
    print "A total of {} combinations".format(total)
    print "============="

def run_test2():
    size = 2
    sops = make_bindings(size)
    print '----normal'
    print_sops(sops)
    sops_c = make_bindings(size, on_channel=True, thresh_adjust=0)
    print '----constitutive'
    print_sops(sops_c)

def run_two_modules():
    double_soplist = {}
    sops = make_bindings(3, on_channel=True, thresh_adjust=0)
    keys = sops.keys()
    # Use the cross-product to generate all combinations
    for m1, m2 in itertools.product(keys, keys):
        # Give us DNF form
        both = simplify_logic(m1 | m2, 'dnf')
        # Some of them are repeated
        b = double_soplist.setdefault(both, [])
        # Don't append the keys, append the amounts
        total = len(sops[m1]) * len(sops[m2])
        b.append((total))

    print_sops(double_soplist, count=True)

    # print "There are {} items!".format(len(double_soplist))
    # for sop, pairs in double_soplist.items():
    #     print term_to_simple_text(sop), len(pairs)
    # prin
    # one = sops.keys()[0]
    # two = sops.keys()[1]
    # print one, two
    # print one | two

if __name__ == '__main__':
    # run_test2()
    run_two_modules()




# Mutation idea ---
# Mutations add or subtract strengths. If they go through 0, they can turn
# into another kind of binding.. That is cute.

