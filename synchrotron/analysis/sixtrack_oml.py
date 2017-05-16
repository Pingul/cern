import sys
import lossmap as lm
import numpy as np

def hits_from_FirstImpact(file_path):
    ids, turns, colls = np.loadtxt(file_path, dtype=int, skiprows=1, usecols=(0, 1, 2), unpack=True)
    m = colls == 9 # TCP IR3
    return list(zip(turns[m], ids[m]))


if __name__ == "__main__":
    d = sys.argv[1]
    hits = hits_from_FirstImpact(d + "/FirstImpacts.dat")
    hm = lm.CHitMap.from_hits(hits)
    lm.plot(hm)
