import sys
import lossmap as lm
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def hits_from_FirstImpact(file_path):
    ids, turns, colls = np.loadtxt(file_path, dtype=int, skiprows=1, usecols=(0, 1, 2), unpack=True)
    m = colls == 9 # TCP IR3
    return list(zip(turns[m], ids[m]))

def plot_x(filepath):
    x, px = np.loadtxt(filepath, usecols=(0, 1), unpack=True)

    fig, ax = plt.subplots()
    ax.scatter(x, px)
    ax.set_xlabel("x")
    ax.set_ylabel("x'")
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray')
    ax.xaxis.grid(color='gray')
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.1E}".format(x)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.1E}".format(x)))
    plt.show()


if __name__ == "__main__":
    d = sys.argv[1]
    if sys.argv[2] == "export":
        plot_x("../export.txt")
    else:
        plot_x(d + "/dist0.dat")
    # hits = hits_from_FirstImpact(d + "/FirstImpacts.dat")
    # hm = lm.CHitMap.from_hits(hits)
    # lm.plot(hm)
