import sys, os
sys.path.insert(0, "..")

import sixtrack_oml
import matplotlib.pyplot as plt
import lossmap as lm

b1 = sixtrack_oml.Batch("/Users/swretbor/Workspace/work_afs/sixtrack/simulations/oml_study/v16")
b2 = sixtrack_oml.Batch("/Users/swretbor/Workspace/work_afs/sixtrack/simulations/oml_study/v17")

print("Max beam 1/max beam 2:", b1.hitmap.losses(integrated=True).max()/b2.hitmap.losses(integrated=True).max())
lm.plot([b1.hitmap, b2.hitmap], ["SixTrack, beam 1", "SixTrack, beam 2"])
