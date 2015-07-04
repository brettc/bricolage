from sympy.logic import SOPform, simplify_logic
from sympy import symbols, sympify
from sympy.printing.latex import latex

import itertools

def detex(term):
    rt = latex(term)
    rt = rt.replace(r"\left", "")
    rt = rt.replace(r"\right", "")
    rt = rt.replace(r"\wedge", "&")
    rt = rt.replace(r"\vee", "|")
    rt = rt.replace(r"\neg ", "~")
    rt = " ".join(rt.strip().split())
    return rt

def boolean_func_from_coop_binding(factory, channels, bindings):
    """Convert a coop binding into a boolean function"""
    # Can't assume all sites are unique, so we need to reconstruct a the truth
    # table entries to pool together sites that are the same.
    unique = list(set(channels))
    indexes = {unique[i]: i for i in range(len(unique))}

    # Generate all possible states 
    all_states = [_ for _ in itertools.product([0, 1], repeat=len(bindings))]

    # We'll put the truth table entries that are TRUE into this
    is_true = []
    for state in all_states:
        # We squash our entries down to the unique entries
        entry = [0] * len(unique)
        summed = 0
        for s, c, b in zip(state, channels, bindings):
            # If the current column entry is true...
            if s:
                # ... make sure the corresponding entry is true
                entry[indexes[c]] = 1
                # And add the binding strength in
                summed += b

        # This tracks what we do in the cpp class ...
        if summed >= len(bindings):
            is_true.append(entry)

    # This should be possible, but sympy barfs on 'E1' for some bizarre reason
    # names = [factory.name_for_channel(u) for u in unique] 
    # So, we use simple names and the replace at the end...
    names = list('abcdefghijklmnopqrstuvwxyz')[:len(unique)]
    sop = SOPform(names, is_true)

    # Now force the 0 and 1 (on and off channels) to be evaluated
    for i, u in enumerate(unique):
        # off channel
        if u == 0:
            # This is ALWAYS OFF
            sop = sop.subs(names[i], 0)
        # On channel
        if u == 1:
            # This is ALWAYS ON
            sop = sop.subs(names[i], 1)

    # Simplify the logic again
    sop = simplify_logic(sop)

    if sop == False:
        return "OFF"
    if sop == True:
        return "ON"

    text = detex(sop)

    # Necessary cos of the E1 problem
    for i, n in enumerate(names): 
        text = text.replace(n, factory.name_for_channel(unique[i]))

    return text

