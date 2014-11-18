"""
A Quick test of how to convert cooperative binding strengths into boolean
operators. This generates (roughly) the kind of boolean operators that we see
in the Davidson-Peters models. IE. No crazy NAND, XOR or anything.
"""


import itertools
from sympy.logic import SOPform, POSform
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

    strengths = range(-size, size+1)
    letters = [chr(ord('a') + i) for i in range(size)]

    soplist = {}
    threshold = size - thresh_adjust

    inputs = list(itertools.product((0, 1), repeat=size))
    for bindings in itertools.product(strengths, repeat=size):
        minterms = []
        for i in inputs:
            active = [a for a, b in zip(bindings, i) if b == 1]
            if sum(active) >= threshold:
                minterms.append(i)
        sop = the_form(letters, minterms)
        if len(sop.atoms()) > 0:
            b = soplist.setdefault(sop, [])
            b.append(bindings)

    # for sop, bind in soplist.items():
    #     print detex(sop), len(bind)
    #     # for b in bind:
    #     #     print b,
    #     # print
            
    return soplist

# sops = make_bindings(3, thresh_adjust=-1)
sops = make_bindings(3)

# Look at the distribution of boolean functions that are generated
sort_sops = [(len(v), k) for k, v in sops.items()]
for p in sorted(sort_sops, key=lambda x: x[0]):
    print p


