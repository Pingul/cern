
import sys, os
sys.path.insert(0, "..")

from lhccomp import LHCComparison
from phasespace import PhaseSpace
import lossmap as lm
import oml
import analyse_fill as af

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.optimize import least_squares

from settings import settings


dir_in = "/Users/swretbor/Workspace/collimation/proj/synchrotron/simulations/store/pdf_inside"
dir_out = "/Users/swretbor/Workspace/collimation/proj/synchrotron/simulations/store/pdf_outside"
# dir_in = "/Users/swretbor/Workspace/collimation/proj/synchrotron/simulations/store/pdf_inside_1turnmap"
# dir_out = "/Users/swretbor/Workspace/collimation/proj/synchrotron/simulations/store/pdf_outside_1turnmap"

fill = af.aggregate_fill(1, from_cache=True)

inside = {
        "ps" : PhaseSpace(dir_in + "/" + settings.STARTDIST_FILE),
        "lm" : lm.CHitMap(dir_in + "/" + settings.COLL_FILE)
        }
ic = LHCComparison(fill, inside['ps'], inside['lm'])

outside = {
        "ps" : PhaseSpace(dir_out + "/" + settings.STARTDIST_FILE),
        "lm" : lm.CHitMap(dir_out + "/" + settings.COLL_FILE)
        }
oc = LHCComparison(fill, outside['ps'], outside['lm'])

# ic.halign()
# oc.halign()

tmin = min(ic.t().min(), oc.t().min())
tmax = max(ic.t().max(), oc.t().max())

t = np.arange(tmin, tmax, 0.1)
y1 = interpolate.interp1d(ic.t(), ic.BLM(False), fill_value=0.0, bounds_error=False)(t)
y2 = interpolate.interp1d(oc.t(), oc.BLM(False), fill_value=0.0, bounds_error=False)(t)
fy = interpolate.interp1d(*fill.blm_ir3())(t)
y = np.stack((y1, y2))

model = lambda c: np.sum(c.reshape(2, 1)*y, axis=0)
error = lambda c: np.sum((model(c) - fy)**2)
c0 = np.array([1.0, 1.0])

res = least_squares(error, c0, loss='cauchy', bounds=([0.0, np.inf]), max_nfev=50000, jac='3-point')
print(res)
print("Inside: ", res.x[0]/sum(res.x))
print("Outside: ", res.x[1]/sum(res.x))

fig, ax = plt.subplots()
ax.plot(*fill.blm_ir3(), label="Aggregate fill", color='r')
ax.plot(ic.t(), ic.BLM(False), label='inside', color='orange')
ax.plot(oc.t(), oc.BLM(False), label='outside', color='purple')
ax.plot(t, model(c0), label='combined', zorder=5, linestyle='--', color='b')
ax.plot(t, model(res.x), label='optimal', linestyle='--', color='g')
ax.legend(loc="upper left")
ax.set_xlabel("t (s)")
ax.set_ylabel("BLM signal")
ax.set_yscale("log")
ax.set_xlim([0, fill.crossover_point()['t']])
plt.show()
