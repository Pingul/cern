import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

import numpy as np
import sys

DIR = sys.argv[1]
DYNKSETS = DIR + '/dynksets.dat'
PARTICLE_ENERGY = "part_1_energy.txt"

dynk_val = []

with open(DYNKSETS, 'r') as f:
    for line in f.readlines():
        if line.lstrip().startswith("#") or not line: continue
        turn, element, attrib, setid, funname, value = line.rstrip().split()
        dynk_val.append(float(value))

dynk_val = np.array(dynk_val)
size = len(dynk_val)

total_energy = np.empty(size)
total_momentum = np.empty(size)

with open(PARTICLE_ENERGY, 'r') as f:
    for i, line in enumerate(f.readlines()):
        momentum, tot_e = map(float, line.strip().split())
        total_energy[i] = tot_e
        total_momentum[i] = momentum

turns = range(size)

fig, ax = plt.subplots()
ax.plot(turns, dynk_val, color='black', label='ref')
ax.plot(turns, total_energy, color='r', label='total energy')
ax.plot(turns, total_momentum, color='b', label='momentum')
ax.legend(loc="upper left")

ax.set_xlabel("Turns")
ax.set_ylabel("Energy (GeV)")
ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e3)))
fig.suptitle("Compare Sixtrack momentum and energy")

plt.show()