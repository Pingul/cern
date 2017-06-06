import sys, os
sys.path.insert(0, "..")

import matplotlib.pyplot as plt
from phasespace import PhaseSpace
import lhccomp as lc
import lossmap as lm
import analyse_fill as af
from settings import settings

dirs = [
    # "/Users/swretbor/Workspace/work_afs/2dsynch/slope_test/CInside_C",
    "/Users/swretbor/Workspace/work_afs/2dsynch/slope_test/CInside_C_beta",
    # "/Users/swretbor/Workspace/work_afs/2dsynch/slope_test/CInside_LD",
    "/Users/swretbor/Workspace/work_afs/2dsynch/slope_test/CInside_LD_beta",
    # "/Users/swretbor/Workspace/work_afs/2dsynch/slope_test/CInside_LI",
    "/Users/swretbor/Workspace/work_afs/2dsynch/slope_test/CInside_LI_beta",
]

labels = [
    # "CInside_C", 
    "Constant", 
    # "CInside_LD", 
    "Lin. decrease", 
    # "CInside_LI", 
    "Lin. increase"
]

pss = []
hms = []
for d in dirs:
    print("Reading '{}'".format(d))
    hm = lm.CHitMap(d + "/" + settings.COLL_FILE)
    hms.append(hm)
    pss.append(PhaseSpace(d + "/" + settings.STARTDIST_FILE))


print("Plot")

fill = af.aggregate_fill(1, from_cache=True)

fig, loss_ax = plt.subplots()
loss_ax.plot(*fill.blm_ir3(), label="Aggr. fill (beam {})".format(fill.beam), color='r')
loss_ax.set_yscale("log")
loss_ax.set_xlim([-5, 40])
loss_ax.set_ylim([0.5e-5, 1])
loss_ax.set_ylabel("Losses (∆particles/1.3s)")
loss_ax.set_xlabel("t (s)")
loss_ax.axvspan(0.0, fill.crossover_point()['t'], facecolor='b', zorder=0, alpha=0.1)

for i, hm in enumerate(hms):
    comp = lc.LHCComparison(fill, None, hm)
    loss_ax.plot(comp.t(), comp.BLM(), label=labels[i])

loss_ax.legend(loc="upper right")

plt.title("Compare 2d-synchrotron with aggregate fill")
plt.draw()

fig, ax = plt.subplots()
s = [(ps.h - settings.H_SEPARATRIX).astype(int) for ps in pss]
ax.hist(s, edgecolor='white')
ax.set_title("Action value H starting distribution")
ax.set_xlabel("∆H")
ax.set_ylabel("#")
plt.show()
