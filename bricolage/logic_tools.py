from sympy.logic import SOPform, simplify_logic
from sympy import symbols
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
        for i, s in enumerate(state):
            # Ignore states where CHANNEL 0 is on
            if s and channels[i] == 0:
                continue
            # or CHANNEL 1 is off
            if not s and channels[i] == 1:
                continue

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

        if summed >= len(bindings):
            is_true.append(entry)

    # For some bizarre I don't yet understand, the variable names sometimes
    # cause it to crash
    names = list('abcdefghijklmnopqrstuvwxyz')[:len(unique)]

    sop = SOPform(names, is_true)

    # # Can we simplify this, once we know the off/on channels?
    # off_c, on_c = symbols("a,b")
    # off_c = False
    # on_c = True
    # sop = simplify_logic(sop & off_c, form='dnf', deep=True)

    if sop == False:
        return "OFF"
    if sop == True:
        return "ON"

    text = detex(sop)
    for i, n in enumerate(names): 
        text = text.replace(n, factory.name_for_channel(unique[i]))

    return text

