import sys
sys.path.append("../comp_sixtrack_toymodel")
from comp import phi_to_z

import numpy as np
import pandas as pd
import plot as lm
from phasespace import PhaseSpace

SIXTRACK_EXPORT_FILE = "sixtrack.dat"

# Don't really know how these work. Should maybe be reviewed
def export_fitted_ramp():
    ramp = read_ramp(RAMP_FILE, 500000)
    with open("ramp_fitted.txt", 'w') as f:
        for i, v in enumerate(ramp['e_fitted']):
            f.write("{} {:.16f}\n".format(i, v/1e6))

def export_particles(phasespace, plist, save_file):
    print("writing {} particles to '{}'".format(len(plist), save_file))
    with open(save_file, 'w') as f:
        f.write("{},1\n".format(len(plist)))
        for pid in plist:
            f.write("{},{},{}\n".format(phasespace.denergy[pid], phasespace.phase[pid], phasespace.h[pid]))

def distribution_to_sixtrack(dist_path, save_file):
    ps = PhaseSpace(dist_path)
    with open(save_file, 'w') as f:
        for i in range(ps.nbr_p):
            x, px, y, py = 0, 0, 0, 0
            z = phi_to_z(ps.phase[i])
            E = 450e3 + ps.denergy[i]*1e-6 # MeV
            f.write(" {:>17.12f}, {:>17.12f}, {:>17.12f}, {:>17.12f}, {:>17.12f}, {:>17.12f}\n".format(x, px, y, py, z, E))


ACTION = sys.argv[1]
if ACTION == 'sixtrack':
    if not len(sys.argv) > 2: raise Exception("need to provide a path to a distribution to export")
    print("export '{}' to sixtrack".format(sys.argv[2]))
    print("saving data to '{}'".format(SIXTRACK_EXPORT_FILE))
    distribution_to_sixtrack(sys.argv[2], SIXTRACK_EXPORT_FILE)
    print("export finished")
else:
    print("unrecognised action")

