import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from scipy.ndimage.interpolation import shift

import numpy as np
import sys

DIR = sys.argv[1]
DUMP = DIR + '/DUMP.txt'
DYNKSETS = DIR + '/dynksets.dat'
P1 = 'p1.txt'
P2 = 'p2.txt'

dynk_val = []

with open(DYNKSETS, 'r') as f:
    for line in f.readlines():
        if line.lstrip().startswith("#") or not line: continue
        turn, element, attrib, setid, funname, value = line.rstrip().split()
        dynk_val.append(float(value))

dynk_val = np.array(dynk_val)
size = len(dynk_val)

p1 = {'e_real' : np.empty(size), 'e_delta' : np.empty(size)}
p2 = {'e_real' : np.empty(size), 'e_delta' : np.empty(size)}

with open(P1, 'r') as f:
    for i, line in enumerate(f.readlines()):
        id, e = map(float, line.rstrip().split())
        p1['e_real'][i] = e

with open(P2, 'r') as f:
    for i, line in enumerate(f.readlines()):
        id, e = map(float, line.rstrip().split())
        p2['e_real'][i] = e

with open(DUMP, 'r') as f:
    it1 = it2 = 0
    for line in f.readlines():
        if line.lstrip().startswith("#"): continue

        id, turn, s, x, xp, y, yp, z, de, ktrack = line.rstrip().split()
        if int(id) == 1:
            p1['e_delta'][it1] = float(de)
            it1 += 1
        else:
            p2['e_delta'][it2] = float(de)
            it2 += 1

def mean_sq_err(x, y):
    if len(x) != len(y):
        print(len(x), len(y))
        raise Exception("x and y needs to have the same dimension")
    return 1.0/len(x) * sum((x - y)**2)

def fit(x, test_func, param=1.0, threshold=0.00001, maxiter=20):
    print("--- FIT ---")
    err = test_func(x)
    n = 0
    step = param
    while err > threshold and n < maxiter:
        x_tmp = x + param
        err_tmp = test_func(x_tmp)
        deriv = err_tmp - err
        if deriv < 0: # we got closer
            param += step
            err = err_tmp
        else:
            param -= step
            step = -step/2.0
            param += step
        print("mean_sq: ", err_tmp, " d: ", deriv, "p: ", param)

        n += 1
    else:
        print("--- default param within error ---")
        return 0.0
    print("--- BEST FIT: {} ---".format(param))
    return param

p1['e_calc'] = dynk_val*(1.0 + p1['e_delta'])
# p1['e_calc'] = shift(p1['e_calc'], 1, cval=450e3)

# fitting energy
# p = fit(p1['e_calc'], lambda x : mean_sq_err(x, p1['e_real']))
# p1['e_calc'] += p

mean_sq = mean_sq_err(p1['e_calc'], p1['e_real'])
print("Error: ", mean_sq)


# Plot
fig, ax = plt.subplots()

turns = range(size)
ax.plot(turns, dynk_val, color='black', zorder=1)
ax.plot(turns, p1['e_real'], color='r', zorder=2)
ax.plot(turns, p1['e_calc'], color='b', linestyle='--', zorder=3)

ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e3)))
ax.set_ylabel("Energy (GeV)")
ax.set_xlabel("Turns")

plt.show()