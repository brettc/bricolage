from .logic2 import Parameters
from .threshold3_ext import Factory

p = Parameters()
f = Factory(p)
n = f.create_network()

for g in n.genes:
    print g.modules



