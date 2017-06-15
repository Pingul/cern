import sys, os
import lossmap as lm
import phasespace as ps
import lhccomp
import plot as psplot
import numpy as np
import pickle
import warnings
import analyse_fill as af

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from settings import settings
from logger import ModuleLogger, LogLevel
lg = ModuleLogger("sixbatch")

class Batch:
    def __init__(self, path, forced=False):
        self.path = path
        self.nbr_jobs = 0
        self.nbr_valid_jobs = 0

        self.phasespace = None
        self.hits = []
        self.hitmap = None

        if path:
            cache_path = "{}/{}".format(path, settings.CACHE_BATCH_FILE)
            if not forced and os.path.exists(cache_path):
                lg.log("trying to read from cache")
                self.load_cache(cache_path)
            else:
                lg.log("cache file could not be found or forced to fetch", log_level=LogLevel.notify)
                self.aggregate()
                self.cache(cache_path)
        else:
            self.aggregate()
    
    def pack(self):
        return {
            'path' : self.path,
            'nbr_jobs' : self.nbr_jobs,
            'nbr_valid_jobs' : self.nbr_valid_jobs,
            'phasespace' : self.phasespace,
            'hitmap' : self.hitmap
        }

    def unpack(self, dump):
        self.path = dump['path']
        self.nbr_jobs = dump['nbr_jobs']
        self.nbr_valid_jobs = dump['nbr_valid_jobs']
        self.phasespace = dump['phasespace']
        self.hitmap = dump['hitmap']

    def cache(self, file_path):
        lg.log("caching batch at '{}'".format(file_path))
        with open(file_path, 'wb') as f:
            pickle.dump(self.pack(), f)

    def load_cache(self, file_path):
        lg.log("loading batch from '{}'".format(file_path))
        with open(file_path, 'rb') as f:
            self.unpack(pickle.loads(f.read()))

    def aggregate(self):
        n = 1
        required_files = ["chits.txt", "DUMP.txt", "fort.13"]
        lg.log("aggregating batch from '{}'".format(self.path))
        print("{}/{}{:04d}".format(self.path, "run", n))
        while os.path.isdir("{}/{}{:04d}".format(self.path, "run", n)):
            job_path = "{}/{}{:04d}".format(self.path, "run", n)
            lg.log("trying '{}'...".format(job_path), end=" ")
            self.nbr_jobs += 1
            for req in required_files:
                if not os.path.exists("{}/{}".format(job_path, req)):
                    lg.log("not ok", module_prestring=False, log_level=LogLevel.warning)
                    break;
            else:
                # pid_offset = self.phasespace.nbr_p if self.phasespace else 0
                pid_offset = self.hitmap.ids.size if self.hitmap else 0
                phasespace = ps.PhaseSpace.from_fort13("{}/{}".format(job_path, "fort.13"))
                # hm = hits_from_FirstImpact("{}/{}".format(job_path, "FirstImpacts.dat"), pid_offset)
                hm = hits_from_chits("{}/{}".format(job_path, "chits.txt"), pid_offset)
                # hm = lm.CHitMap("{}/{}".format(job_path, settings.COLL_FILE), pid_offset=pid_offset, store_ip=False)
                if self.phasespace is None:
                    self.phasespace = phasespace
                    self.hitmap = hm
                else:
                    self.phasespace = ps.PhaseSpace.merge_two(self.phasespace, phasespace)
                    self.hitmap = self.hitmap.concatenate(hm)
                self.nbr_valid_jobs += 1
                lg.log("ok", module_prestring=False, log_level=LogLevel.success)
            n += 1
            

def hits_from_FirstImpact(filepath, offset):
    try:
        ids, turns, colls = np.loadtxt(filepath, dtype=int, skiprows=1, usecols=(0, 1, 2), unpack=True)
    except ValueError:
        return lm.CHitMap(store_ip=False)
    m = colls == 9 # TCP IR3
    return lm.CHitMap.from_hits(list(zip(turns[m], ids[m])), pid_offset=offset)

def hits_from_chits(filepath, offset):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # will generate warnings if no hits
        turns = np.loadtxt(filepath, dtype=int, skiprows=1, usecols=(1,), unpack=True).reshape(-1) # In case of only 1 hit, this becomes a 0-d scalar and not a proper array https://stackoverflow.com/questions/773030/why-are-0d-arrays-in-numpy-not-considered-scalar. Reshaping solves the issue.
        coll = np.loadtxt(filepath, dtype=np.object, skiprows=1, usecols=(0,), unpack=True).reshape(-1)

    # print([i for i, item in enumerate(coll) if "tcp.6" in item])
    # if np.any(coll == "b'tcp.c6l7.b1'"): print("has loss")

    # creating an ID array here. Current output does not contain ID data,
    # so this is gibberish for now
    ids = np.arange(0, turns.size).astype(int)
    return lm.CHitMap.from_hits(list(zip(turns, ids)), pid_offset=offset)

if __name__ == "__main__":
    action = sys.argv[1]
    if action == "batch":
        d = sys.argv[2]
        b = Batch(d)
        lg.log("Batch stats\n\tJobs: {}\n\tParticles: {}".format(b.nbr_jobs, b.phasespace.nbr_p))
        
        turn_max = b.hitmap.losses(integrated=True).argmax()
        lg.log("Peak at {} s".format(turn_max/11200.0))

        fill = af.aggregate_fill(1, from_cache=True)
        comp = lhccomp.LHCComparison(fill, b.phasespace, b.hitmap)
        lhccomp.plot_comp(fill, (comp.t(), comp.BLM()), block=False)
        b.phasespace.plot_particles()
    elif action == "dump":
        phasespace = ps.PhaseSpace.from_six_dump(sys.argv[2])
        lg.log("φ-range: [{:.3f}, {:.3f}]".format(phasespace.phase.min(), phasespace.phase.max()))
        lg.log("E-range: [{:.3E}, {:.3E}]".format(phasespace.denergy.min(), phasespace.denergy.max()))
        phasespace.plot_particles()
    elif action == "dumpx":
        phasespace = ps.PhaseSpace.from_six_dump(sys.argv[2])
        phasespace.x *= 1e3
        phasespace.px *= 1e3
        lg.log(" x-range: [{:.3f}, {:.3f}] mm".format(phasespace.x.min(), phasespace.x.max()))
        lg.log("px-range: [{:.3f}, {:.3f}] mm".format(phasespace.px.min(), phasespace.px.max()))
        ps.plot_x(phasespace)
    else:
        lg.log("unrecognised action")
