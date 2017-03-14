
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
from scipy.optimize import nnls 
from sklearn import linear_model
from bisect import bisect_left

BLM_INT = int(1.3*11245.0)
H_SEPARATRIX = 1909859.317103239
# BLM_INT = int(1)

def compare_to_aggregate():
    tm_lossmap = lm.get_lossmap(lm.COLL_FILE)

    turns = range(min(tm_lossmap.keys()) - BLM_INT, max(tm_lossmap.keys()) + BLM_INT)
    secs = np.array([turn/11245.0 for turn in turns])
    # ramp = np.array(lm.read_ramp(lm.RAMP_FILE, len(turns))['e'])/1.0e9

    losses = np.array([len(tm_lossmap[turn]) if turn in tm_lossmap else 0 for turn in turns])
    synch_avg_losses = lm.moving_average(losses, BLM_INT)

    aggr_fill = oml.aggregate_fill(from_cache=True)
    secs = align(secs, synch_avg_losses, aggr_fill)


    # fitting y = coef*x
    ps = lm.PhaseSpace(lm.STARTDIST_FILE)
    lossmaps, a_values = lm.separate_lossmap(tm_lossmap, ps)
    a_values = np.array(a_values)
    a_values -= round(H_SEPARATRIX)

    prune = True
    if prune:
        H_threshold = -8000
        to_delete = np.where(a_values < H_threshold)
        lossmaps = np.delete(lossmaps, to_delete)
        a_values = np.delete(a_values, to_delete)

    x = []
    for lossmap in lossmaps:
        losses = np.array([len(lossmap[turn]) if turn in lossmap else 0 for turn in turns])
        # avg_loss = lm.moving_average(losses, BLM_INT)
        avg_loss = losses
        x.append(avg_loss)

    # Note: we fit in log space
    xt = np.transpose(x)
    y = interpolate.interp1d(*aggr_fill.oml())(secs)

    method = 4
    if method == 1:
        print("log:sklearn")
        xt[xt<=0] = sys.float_info.min
        xt = np.log(xt)
        y = np.log(y)

        clf = linear_model.LinearRegression()
        clf.fit(xt, y)
        coef = clf.coef_
    elif method == 2:
        print("log:nnls")
        xt[xt<=0] = sys.float_info.min
        xt = np.log(xt)
        y = np.log(y)

        coef = nnls(xt, y)[0]
    elif method == 3:
        print("log:numpy")
        xt[xt<=0] = sys.float_info.min
        xt = np.log(xt)
        y = np.log(y)

        coef = np.linalg.lstsq(xt, y)[0]
    elif method == 4:
        print("nnls")
        coef = nnls(xt, y)[0]

    print("coefficients")
    for i, c in enumerate(coef):
        print("\t{:>5} = {:<10.3f} (H = {})".format("a{}".format(i), c, a_values[i]))
    x_fit = np.sum(coef*xt, axis=1)
    x_fit = lm.moving_average(x_fit, BLM_INT)

    if method <= 3:
        x_fit = np.exp(x_fit)


    # Plotting lossmap fit
    fig, loss_ax = plt.subplots()

    loss_ax.plot(secs, x_fit,  color='green', label="least square fit", linestyle='--', zorder=6)
    loss_ax.plot(secs, synch_avg_losses, color='orange', zorder=5, label='2d-synchrotron')
    loss_ax.plot(*aggr_fill.oml(), color='red', label='LHC data')
    loss_ax.set_ylabel("Losses (∆particles/1.3s)")
    loss_ax.set_xlabel("t (s)")
    loss_ax.set_yscale('log')
    loss_ax.axvspan(0.0, 13.55, facecolor='b', zorder=0, alpha=0.1)
    loss_ax.legend(loc="upper right")

    e_ax = loss_ax.twinx()
    e_ax.set_xlabel("t (s)")
    e_ax.plot(*aggr_fill.data["energy"], color='black')
    e_ax.set_ylabel("E (GeV)")

    loss_ax.set_xlim([-10, 40])
    fig.suptitle("Aggregate vs 2d-synchrotron")

    # Plotting coefficients
    fig, ax = plt.subplots()
    ax.bar(range(len(coef)), coef)
    ax.set_xticks(range(len(coef)))
    ax.set_xticklabels(a_values, rotation='vertical')
    ax.set_ylabel("Value")
    ax.set_xlabel("∆H")
    index = bisect_left(a_values, 0)
    ax.axvspan(index, index + 0.05, facecolor='r', zorder=0, alpha=0.2)
    ax.text(index, max(coef)/2.0, "Separatrix", fontsize=10, 
            ha='center', va='center', rotation='vertical', color='r')
    
    plt.title("Optimized action value coefficients")
    plt.tight_layout()

    plt.show()

def align(secs, losses, aggr_fill):
    """ Aligns losses to aggr_fill by returning a new 'secs' array """

    vloss_peak, iloss_peak = oml.imax(losses)
    vfill_peak, ifill_peak = oml.imax(aggr_fill.oml().y)
    delta = secs[iloss_peak] - aggr_fill.oml().x[ifill_peak]
    print("peaks\n\taggregate: {:.2f}\n\t2d-synch : {:.2f}\n\tdelta    : {:.2f}"
            .format(secs[iloss_peak], aggr_fill.oml().x[ifill_peak], delta))
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
