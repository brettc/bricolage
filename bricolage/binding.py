"""
A Quick test of how to convert cooperative binding strengths into boolean
operators. This generates (roughly) the kind of boolean operators that we see
in the Davidson-Peters models. IE. No crazy NAND, XOR or anything.
"""


import itertools
from sympy.logic import SOPform #, POSform
the_form = SOPform
from sympy.printing.latex import latex

def detex(term):
    rt = latex(term)
    rt = rt.replace(r"\left", "")
    rt = rt.replace(r"\right", "")
    rt = rt.replace(r"\wedge", "&")
    rt = rt.replace(r"\vee", "|")
    rt = rt.replace(r"\neg ", "~")
    rt = rt.replace(" ", "")
    return rt

def make_bindings(size, thresh_adjust=0):
    """
    The rough idea. Each "site" either contributes to or detracts from the
    binding. Assuming there are three sites, We assign each site from -3 to +3.
    If a site is "bound", then we add the amount at that site. 
    """

    strengths = range(-size, size+1)
    # Assign each site a..z
    letters = [chr(ord('a') + i) for i in range(size)]

    soplist = {}
    threshold = size - thresh_adjust

    # This generates all possible combinations of inputs for the number binary
    # sites
    inputs = list(itertools.product((0, 1), repeat=size))

    # This generates all possible combinations of binding strengths.
    for bindings in itertools.product(strengths, repeat=size):
        minterms = []

        # Iterate through all inputs for this
        for i in inputs:
            active = [a for a, b in zip(bindings, i) if b == 1]
            if sum(active) >= threshold:
                minterms.append(i)

        # Now, what boolean function does that set of combinations work for?
        sop = the_form(letters, minterms)

        if len(sop.atoms()) > 0:
            # Some of them are repeated
            b = soplist.setdefault(sop, [])
            b.append(bindings)

    return soplist

sops = make_bindings(3)

# Look at the distribution of boolean functions that are generated. Some occur
# more than others ... We might want to change these via probabilities..?
sort_sops = [(len(v), k) for k, v in sops.items()]
for p in sorted(sort_sops, key=lambda x: x[0]):
    print p


# Mutation idea ---
#
# Mutations add or subtract strengths. If they go through 0, they can turn
# into another kind of binding.. That is cute.

