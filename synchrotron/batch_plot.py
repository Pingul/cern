import pickle
from phasespace import PhaseSpace
from lhccomp import fit_to_LHC_aggregate, compare_to_LHC_aggregate
from lossmap import *
from settings import *
from plot import plot_hamiltonian_dist_histogram

import os, sys
sys.path.append("../common/")
from logger import ModuleLogger, LogLevel
lg = ModuleLogger("batch")

class Batch:
    def __init__(self, path, forced=settings.BATCH_FORCE_FETCH):
        self.path = path
        self.nbr_jobs = 0
        self.nbr_valid_jobs = 0

        self.ps = None
        self.hits = []
        self.lossmap = None

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
            'lossmap' : self.lossmap
        }

    def unpack(self, dump):
        self.path = dump['path']
        self.nbr_jobs = dump['nbr_jobs']
        self.nbr_valid_jobs = dump['nbr_valid_jobs']
        self.ps = dump['ps']
        self.lossmap = dump['lossmap']

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
        required_files = ["stdout.txt", "startdist.dat", "coll.dat"]
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
                # try:
                pid_offset = self.ps.nbr_p if self.ps else 0
                job_ps = PhaseSpace("{}/startdist.dat".format(job_path), mute=True)
                self.ps = PhaseSpace.merge_two(self.ps, job_ps)
                
                hits = hits_from_collfile("{}/coll.dat".format(job_path), pid_offset=pid_offset, mute=True)
                self.hits += hits
    
                self.nbr_valid_jobs += 1
                lg.log("ok", module_prestring=False, log_level=LogLevel.success)
                # except:
                    # lg.log("not ok")
            n += 1
        self.lossmap = lossmap_from_hits(self.hits)
            

if __name__ == "__main__":
    ACTION = sys.argv[1]
    BATCH_DIR = sys.argv[2] if len(sys.argv) > 2 else settings.BATCH_DIR
    b = Batch(BATCH_DIR)

    if ACTION == "stats":
        lg.log("jobs: {}, {:.1f}% valid".format(b.nbr_jobs, 100*b.nbr_valid_jobs/b.nbr_jobs))
        if b.ps:
            nbr_lost = sum(map(len, b.lossmap.values()))
            lg.log("{} particles, {:.1f}% lost".format(b.ps.nbr_p, 100*nbr_lost/b.ps.nbr_p))
    elif ACTION == "lossmap":
        plot_lossmap([b.lossmap])
    elif ACTION == "separated-lossmap":
        plot_lossmap(*separate_lossmap(b.lossmap, b.ps))
    elif ACTION == "compare":
        compare_to_LHC_aggregate(b.ps, b.lossmap)
    elif ACTION == "fit":
        fit_to_LHC_aggregate(b.ps, b.lossmap)
    elif ACTION == "startdist":
        lg.log("categorizing particles")
        pbin = b.ps.categorize_particles(b.lossmap)
        lg.log("plotting")
        b.ps.plot_particles(pbin)
    elif ACTION == "ham-dist":
        plot_hamiltonian_dist_histogram(b.ps)
