import pickle
import os, sys

import numpy as np
import lossmap as lm
from phasespace import PhaseSpace
from lhccomp import LHCComparison, plot_comp
from settings import settings
from plot import plot_hamiltonian_dist_histogram, plot_e_dist_histogram
from logger import ModuleLogger, LogLevel
import analyse_fill as af
lg = ModuleLogger("batch")

class Batch:
    def __init__(self, path, forced=settings.BATCH_FORCE_FETCH):
        self.path = path
        self.nbr_jobs = 0
        self.nbr_valid_jobs = 0

        self.ps = None
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
            'ps' : self.ps,
            'hitmap' : self.hitmap
        }

    def unpack(self, dump):
        self.path = dump['path']
        self.nbr_jobs = dump['nbr_jobs']
        self.nbr_valid_jobs = dump['nbr_valid_jobs']
        self.ps = dump['ps']
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
        required_files = ["stdout.txt", settings.STARTDIST_FILE, settings.COLL_FILE, "TCPc_IR7.ch"]
        lg.log("aggregating batch from '{}'".format(self.path))
        while os.path.isdir("{}/{}{}".format(self.path, settings.BATCH_JOB_PRESTRING, n)):
            job_path = "{}/{}{}".format(self.path, settings.BATCH_JOB_PRESTRING, n)
            lg.log("trying '{}'...".format(job_path), end=" ")
            self.nbr_jobs += 1
            for req in required_files:
                if not os.path.exists("{}/{}".format(job_path, req)):
                    lg.log("not ok", module_prestring=False, log_level=LogLevel.warning)
                    break;
            else:
                pid_offset = self.ps.nbr_p if self.ps else 0
                ps = PhaseSpace("{}/{}".format(job_path, settings.STARTDIST_FILE), mute=True)
                hm = lm.CHitMap("{}/{}".format(job_path, settings.COLL_FILE), pid_offset=pid_offset, store_ip=False)
                if self.ps is None:
                    self.ps = ps
                    self.hitmap = hm
                else:
                    self.ps = PhaseSpace.merge_two(self.ps, ps)
                    self.hitmap = self.hitmap.concatenate(hm)
                self.nbr_valid_jobs += 1
                lg.log("ok", module_prestring=False, log_level=LogLevel.success)
            n += 1
            

if __name__ == "__main__":
    ACTION = sys.argv[1]
    BATCH_DIR = sys.argv[2] if len(sys.argv) > 2 else settings.BATCH_DIR
    b = Batch(BATCH_DIR)

    if ACTION == "stats":
        lg.log("jobs: {}, {:.1f}% valid".format(b.nbr_jobs, 100*b.nbr_valid_jobs/b.nbr_jobs))
        if b.ps:
            nbr_lost = b.hitmap.nbr_lost()
            lg.log("{} particles, {:.1f}% lost".format(b.ps.nbr_p, 100*nbr_lost/b.ps.nbr_p))
            lg.log("Particle distribution")
            lg.log("\t{:.2f}% outside above".format(np.sum((b.ps.h > settings.H_SEPARATRIX)*(b.ps.denergy > 0))/b.ps.nbr_p*100.0), module_prestring=False)
            lg.log("\t{:.2f}% inside separatrix".format(np.sum(b.ps.h < settings.H_SEPARATRIX)/b.ps.nbr_p*100.0), module_prestring=False)
            lg.log("\t{:.2f}% outside below".format(np.sum((b.ps.h > settings.H_SEPARATRIX)*(b.ps.denergy < 0))/b.ps.nbr_p*100.0), module_prestring=False)
    elif ACTION == "lossmap":
        lm.plot(b.hitmap)
    elif ACTION == "separated-lossmap":
        lm.plot(*b.hitmap.split(b.ps))
    elif ACTION == "compare":
        fill = af.aggregate_fill(1, from_cache=1)
        comp = LHCComparison(fill, b.ps, b.hitmap)
        # comp.halign()
        plot_comp(fill, (comp.t(), comp.BLM()))
    elif ACTION == "startdist":
        lg.log("categorizing particles")
        pbin = b.ps.categorize_particles(b.hitmap.old_lossmap())
        lg.log("plotting")
        b.ps.plot_particles(pbin)
    elif ACTION == "ham-dist":
        plot_hamiltonian_dist_histogram(b.ps)
    elif ACTION == "e-dist":
        plot_e_dist_histogram(b.ps)
    else:
        lg.log("unrecognised action", log_level=LogLevel.warning)
