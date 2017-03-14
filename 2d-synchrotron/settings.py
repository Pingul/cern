from math import pi

# DIR = "calc/"
# DIR = "simulations/cache/2017.03.13.14.26.09/"
DIR = "simulations/store/5k_H8000/"

PARTICLE_FILE = DIR + "particles.dat"
LINE_FILE = DIR + "lines.dat"
COLL_FILE = DIR + "coll.dat"
STARTDIST_FILE = DIR + "startdist.dat"
ENDDIST_FILE = DIR + "enddist.dat"
RAMP_FILE = "resources/LHC_ramp.dat"

## 1 bucket
PLOT_FRAME = {
    'x' : [-0.1*pi, 2.1*pi],
    'y' : [-0.4e9, 0.4e9]
}

## 3 buckets
# PLOT_FRAME = {
    # 'x' : [-2*pi, 4*pi],
    # 'y' : [-2e9, 2e9]
# }

## Constants
H_SEPARATRIX = 1909859.317103239 # Taken from 2d-synchrotron.hh
