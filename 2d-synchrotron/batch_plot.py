import pickle
from phasespace import PhaseSpace
from lhccomp import compare_to_LHC_aggregate
from lossmap import *
from settings import *

import os
from sys import argv

class Batch:
    def __init__(self, path, forced=settings.BATCH_FORCE_FETCH):
        self.path = path
        self.nbr_jobs = 0
        self.nbr_valid_jobs = 0

        self.ps = None
        self.hits = []
        self.lossmap = None

        if settings.CACHE_BATCH:
            cache_path = "{}/{}".format(settings.BATCH_DIR, settings.CACHE_BATCH_FILE)
            if not forced and os.path.exists(cache_path):
                print("trying to read from cache")
                self.load_cache(cache_path)
            else:
                print("cache file could not be found or forced to fetch")
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
        print("caching batch at '{}'".format(file_path))
        with open(file_path, 'wb') as f:
            pickle.dump(self.pack(), f)

    def load_cache(self, file_path):
        print("loading batch from '{}'".format(file_path))
        with open(file_path, 'rb') as f:
            self.unpack(pickle.loads(f.read()))

    def aggregate(self):
        n = 1
        required_files = ["stdout.txt", "startdist.dat", "coll.dat"]
        print("aggregating batch from '{}'".format(self.path))
        while os.path.isdir("{}/{}{}".format(self.path, settings.BATCH_JOB_PRESTRING, n)):
            job_path = "{}/{}{}".format(self.path, settings.BATCH_JOB_PRESTRING, n)
            print("trying '{}'...".format(job_path), end=" ")
            self.nbr_jobs += 1
            for req in required_files:
                if not os.path.exists("{}/{}".format(job_path, req)):
                    print("not ok")
                    break;
            else:
                # try:
                pid_offset = self.ps.nbr_p if self.ps else 0
                job_ps = PhaseSpace("{}/startdist.dat".format(job_path), mute=True)
                self.ps = PhaseSpace.merge_two(self.ps, job_ps)
                
                hits = hits_from_collfile("{}/coll.dat".format(job_path), pid_offset=pid_offset, mute=True)
                self.hits += hits
    
                self.nbr_valid_jobs += 1
                print("ok")
                # except:
                    # print("not ok")
            n += 1
        self.lossmap = lossmap_from_hits(self.hits)
            

if __name__ == "__main__":
    ACTION = argv[1]

    if ACTION == "stats":
        b = Batch(settings.BATCH_DIR)
        print("jobs: {}, {:.1f}% valid".format(b.nbr_jobs, 100*b.nbr_valid_jobs/b.nbr_jobs))
        if b.ps:
            nbr_lost = sum(map(len, b.lossmap.values()))
            print("{} particles, {:.1f}% lost".format(b.ps.nbr_p, 100*nbr_lost/b.ps.nbr_p))
    elif ACTION == "lossmap":
        b = Batch(settings.BATCH_DIR)
        plot_lossmap([b.lossmap])
    elif ACTION == "separated-lossmap":
        b = Batch(settings.BATCH_DIR)
        plot_lossmap(*separate_lossmap(b.lossmap, b.ps))
    elif ACTION == "compare":
        b = Batch(settings.BATCH_DIR)
        compare_to_LHC_aggregate(b.ps, b.lossmap)
    elif ACTION == "startdist":
        b = Batch(settings.BATCH_DIR)
        print("categorizing particles")
        pbin = b.ps.categorize_particles(b.lossmap)
        print("plotting")
        b.ps.plot_particles(pbin)
