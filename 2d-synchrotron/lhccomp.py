
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import numpy as np
import pickle

import sys
sys.path.append("/Users/swretbor/Workspace/collimation/proj/lhcstat")
import plot as lm
import oml

from scipy import interpolate

from bisect import bisect_left

BLM_INT = int(1.3*11245.0)

def compare_to_aggregate():
    lossmap = lm.get_lossmap(lm.COLL_FILE)

    turns = range(max(lossmap.keys()))
    secs = np.array([turn/11245.0 for turn in turns])
    ramp = np.array(lm.read_ramp(lm.RAMP_FILE, len(turns))['e'])/1.0e9

    losses = np.array([len(lossmap[turn]) if turn in lossmap else 0 for turn in turns])
    avg_losses = lm.moving_average(losses, BLM_INT)

    aggr_fill = oml.aggregate_fill(from_cache=True)
    secs_aligned = align(secs, avg_losses, aggr_fill)

    # OPTIMIZE
    with open("opt_losses.dat", 'rb') as f:
        opt_losses = pickle.loads(f.read())

    # opt_losses = optimize(secs_aligned, losses, aggr_fill)
    # with open("opt_losses.dat", 'wb') as f:
        # pickle.dump(opt_losses, f)

    # opt_lossmap = regenerate_lossmap(opt_losses, lossmap)
    # print(opt_lossmap)

    # avg_losses = lm.moving_average(opt_losses, BLM_INT)

    # Plotting
    fig, loss_ax = plt.subplots()

    loss_ax.plot(secs_aligned, avg_losses, color='orange', zorder=5, label='2d-synchrotron')
    loss_ax.plot(*aggr_fill.oml(), color='red', label='LHC data')
    loss_ax.set_ylabel("Losses (∆particles/1.3s)")
    loss_ax.set_xlabel("t (s)")
    loss_ax.set_yscale('log')
    loss_ax.axvspan(0.0, 13.55, facecolor='b', zorder=0, alpha=0.1)
    loss_ax.legend(loc="upper right")

    e_ax = loss_ax.twinx()
    e_ax.set_xlabel("t (s)")
    # e_ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1000.0)))

    e_ax.plot(secs_aligned, ramp, color='gray', linestyle='--', zorder=4, label='LHC energy ramp')[0]
    e_ax.plot(*aggr_fill.data["energy"], color='black')
    e_ax.set_ylabel("E (GeV)")

    loss_ax.set_xlim([-10, 40])

    fig.suptitle("Aggregate vs 2d-synchrotron")
    plt.show()

def align(secs, losses, aggr_fill):
    """ Aligns losses to aggr_fill by returning a new 'secs' array """

    vloss_peak, iloss_peak = oml.imax(losses)
    vfill_peak, ifill_peak = oml.imax(aggr_fill.oml().y)
    delta = secs[iloss_peak] - aggr_fill.oml().x[ifill_peak]
    print("shift by", delta, "seconds")
    return secs - delta

def optimize(secs, losses, aggr_fill):
    """ Lets do a naive implementation:
            1. Remove a random particle
            2. Recalculate the lossmap. If one timestamp has lossmap < lhcdata, put back the particle
            3. If we have failed to remove a particle n number of times, stop. Else, continue with #1.
    """

    # We want aggr_fill to use the 'secs' timescale
    oml_y = interpolate.interp1d(*aggr_fill.oml())(secs)
    avg_loss = lm.moving_average(losses, BLM_INT)
    lm_win = np.sum(avg_loss > oml_y)
    
    failed = 0
    removed = 0
    printed = set()
    while failed < 400:
        if removed % 100 == 0 and removed not in printed: 
            print("succesfully removed:", removed)
            printed.add(removed)

        temp_losses = np.array(losses)
        temp_avg_loss = np.array(avg_loss)
        
        r = 0
        while True:
            r = np.random.randint(0, len(temp_losses))
            if temp_losses[r] > 0:
                temp_losses[r] -= 1
                break

        rstart = int(r - BLM_INT/2)
        rend = int(r + BLM_INT/2)
        temp_avg_loss[rstart:rend] -= 1.0/BLM_INT

        if np.sum(temp_avg_loss > oml_y) < lm_win:
            failed += 1
        else:
            losses = temp_losses
            avg_loss = temp_avg_loss
            removed += 1
            failed = 0

    print("threshold reached")
    print("removed", removed, "particles")
    return losses

def regenerate_lossmap(opt_losses, lossmap):
    """ 'opt_losses' should have been optimized. Reflect the trim in 'opt_losses' to the lossmap.  """

    # opt_losses [losses_turn_0, losses_turn_1, ...]
    # lossmap {turn : [id0, id1], ...}


    turns = range(max(lossmap.keys()))
    losses = np.array([len(lossmap[turn]) if turn in lossmap else 0 for turn in turns])
    loss_delta = losses - opt_losses

    reduced_lossmap = lossmap

    s = 0
    for turn in np.where(loss_delta > 0)[0]:
        remove = set()
        while len(remove) < loss_delta[turn]:
            remove.add(np.random.randint(0, len(lossmap[turn])))
        reduced_lossmap[turn] = np.delete(reduced_lossmap[turn], )
    return reduced_lossmap


compare_to_aggregate()
