from sympy.logic import SOPform, simplify_logic
from sympy.printing.latex import latex
from bricolage.core import InterventionState
import itertools


def function_to_text(eq):

    if eq == False:
        return "OFF"
    elif eq == True:
        return "ON"

    text = latex(eq)
    text = text.replace(r"\left", "")
    text = text.replace(r"\right", "")
    text = text.replace(r"\wedge", "&")
    text = text.replace(r"\vee", "|")
    text = text.replace(r"\neg ", "~")
    text = text.replace(r"_", "")
    text = " ".join(text.strip().split())

    return text


def text_for_gene(gene):
    eq = compute_gene_function(gene)
    return function_to_text(eq)


def text_for_cis_mod(mod):
    eq = compute_cis_mod_function(mod)
    return function_to_text(eq)


def compute_gene_function(gene):
    eq = False
    for mod in gene.modules:
        # Skip anything that is turned OFF
        if InterventionState.INTERVENE_OFF == mod.intervene:
            continue
        eq |= compute_cis_mod_function(mod)

    return simplify_logic(eq)


def compute_cis_mod_function(mod):
    wrld = mod.gene.network.factory.world
    unique = list(set(mod.channels))

    # Remove the fixed channels
    if 0 in unique:
        unique.remove(0)
    if 1 in unique:
        unique.remove(1)

    # Build a list of True stats
    is_true = []

    # What if there's nothing here...?
    if not unique:
        st = wrld.create_state()

        # Set our fixed channel
        st[1] = 1
        act = mod.is_active(st)
        if act:
            return True
        return False

    # Generate all possible states for this many channels
    all_states = [_ for _ in itertools.product([0, 1], repeat=len(unique))]
    for state in all_states:
        channel_state = wrld.create_state()
        channel_state[1] = 1
        for i, channel in enumerate(state):
            if channel:
                channel_state[unique[i]] = 1

        # Now, ASK the module whether this state turns it on...
        act = mod.is_active(channel_state)
        if act:
            is_true.append(state)

    # NOTE: We need to append _, as sympy barfs on 'E1' for
    # some bizarre reason...
    names = ['_' + wrld.name_for_channel(u) for u in unique]
    sop = SOPform(names, is_true)

    return sop


# NOTE: keep this around for reference for now.
def old_boolean_func_from_coop_binding(world, channels, bindings):
    """Convert a coop binding into a boolean function"""
    # Can't assume all sites are unique, so we need to reconstruct a truth
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
    # names = [world.name_for_channel(u) for u in unique]
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
        text = text.replace(n, world.name_for_channel(unique[i]))

    return text
