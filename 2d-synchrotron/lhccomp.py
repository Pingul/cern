
import matplotlib.pyplot as plt
import numpy as np

import lossmap as lm
from phasespace import PhaseSpace

import sys, os
sys.path.append("/Users/swretbor/Workspace/collimation/proj/lhcstat")
import oml
import analyse_fill as af

from scipy import interpolate
from scipy.optimize import nnls 
from sklearn import linear_model
from bisect import bisect_left

from settings import *

def trailing_integration(sequence, N):
    """ For each entry, integrate the N first entries.
    """

    sums = np.convolve(sequence, np.ones((N,)))
    return sums[:len(sequence)]

def compare_to_LHC_aggregate(ps, lossmap):
    beam = 1

    turns = np.array(range(max(lossmap.keys())))
    losses = np.array([len(lossmap[turn]) if turn in lossmap else 0 for turn in turns])
    BLM_2dsynch = trailing_integration(losses, settings.BLM_INT)

    aggr_fill = af.aggregate_fill(beam, from_cache=True)
    secs = np.array([turn/11245.0 for turn in turns])
    asecs = halign(secs, BLM_2dsynch, aggr_fill)
    BLM_2dsynch = valign(aggr_fill, asecs, BLM_2dsynch)

    fig, loss_ax = plt.subplots()
    loss_ax.plot(*aggr_fill.blm_ir3(), label="Aggr. fill (beam {})".format(beam), color='r')
    loss_ax.plot(asecs, BLM_2dsynch, label="2d-synch BLM")
    loss_ax.set_yscale("log")
    loss_ax.set_xlim([-5, 40])
    loss_ax.set_ylim([0.5e-5, 1])
    loss_ax.legend(loc="upper right")
    plt.title("Compare 2d-synchrotron with aggregate fill")
    plt.show()

def fit_to_LHC_aggregate(ps, tm_lossmap):
    """ ps : PhaseSpace of starting distribution
        tm_lossmap : lossmap of coll.dat
    """

    beam = 1
    prune = settings.PRUNE_TIMESCALE
    integration = settings.TRAILING_INTEGRATION # integrate the separated lossmaps before optimisation

    if prune:
        print("pruning time scale")
        turns = range(min(tm_lossmap.keys()) - settings.BLM_INT, int(16.5*11245.0))#max(tm_lossmap.keys()) + BLM_INT)
    else:
        print("using full time scale")
        turns = range(min(tm_lossmap.keys()) - settings.BLM_INT, max(tm_lossmap.keys()) + settings.BLM_INT)

    secs = np.array([turn/11245.0 for turn in turns])
    losses = np.array([len(tm_lossmap[turn]) if turn in tm_lossmap else 0 for turn in turns])
    BLM_2dsynch = trailing_integration(losses, settings.BLM_INT)

    aggr_fill = af.aggregate_fill(beam, from_cache=True)
    asecs = halign(secs, BLM_2dsynch, aggr_fill)
    BLM_2dsynch = valign(aggr_fill, asecs, BLM_2dsynch)

    lossmaps, action_values = lm.separate_lossmap(tm_lossmap, ps)
    action_values = np.array(action_values)
    action_values -= round(settings.H_SEPARATRIX)

    x = []
    for lossmap in lossmaps:
        losses = np.array([len(lossmap[turn]) if turn in lossmap else 0 for turn in turns])
        if integration:
            avg_loss = trailing_integration(losses, settings.BLM_INT)
        else:
            avg_loss = losses.astype(float)
        x.append(avg_loss)

    xt = np.transpose(x)
    y = interpolate.interp1d(*aggr_fill.blm_ir3())(asecs)
    coef = nnls(xt, y)[0]
    print("coefficients")
    for i, c in enumerate(coef):
        print("\t{:>5} = {:<10.5f} (H = {})".format("a{}".format(i), c, action_values[i]))
    x_fit = np.sum(coef*xt, axis=1)

    if not integration:
        x_fit = trailing_integration(x_fit, settings.BLM_INT)
        x_fit = x_fit.reshape(len(x_fit), 1)
        correction = nnls(x_fit, y)[0]
        x_fit *= correction

    # Plotting lossmap fit
    option_string = "(integ.)" if integration else ""
    fig, loss_ax = plt.subplots()

    loss_ax.plot(*aggr_fill.blm_ir3(), color='r', label='Aggr. fill (beam {})'.format(beam))
    loss_ax.plot(secs, BLM_2dsynch, zorder=5, label='2d-synch BLM')
    loss_ax.plot(asecs, x_fit, label="least square fit", linestyle='--', zorder=6, color="forestgreen")
    loss_ax.axvspan(0.0, aggr_fill.crossover_point()['t'], facecolor='b', zorder=0, alpha=0.1)
    loss_ax.set_ylabel("Losses (∆particles/1.3s)")
    loss_ax.set_xlabel("t (s)")
    loss_ax.set_yscale('log')
    loss_ax.set_xlim([-5, 40])
    loss_ax.set_ylim([0.5e-5, 1])
    loss_ax.legend(loc="upper right")
    plt.title("Aggregate vs 2d-synchrotron {}".format(option_string))

    # Plotting coefficients
    fig, ax = plt.subplots()
    ax.bar(range(len(coef)), coef)
    ax.set_xticks(range(len(coef)))
    ax.set_xticklabels(action_values, rotation='vertical')
    ax.set_ylabel("Value")
    ax.set_xlabel("∆H")
    index = bisect_left(action_values, 0)
    ax.axvspan(index, index + 0.05, facecolor='r', zorder=0, alpha=0.2)
    ax.text(index, max(coef)/2.0, "Separatrix", fontsize=10, 
            ha='center', va='center', rotation='vertical', color='r')
    
    plt.title("Optimized action value coefficients {}".format(option_string))
    plt.tight_layout()

    plt.show()

def valign(aggregate_fill, aligned_secs, BLM_2dsynch):
    """ Returns a vertically shifted BLM_2dsynch array that is fitted 
        w.r.t the aggregate fill. We assume that they are horizontally aligned
    """
    BLM_2dsynch = BLM_2dsynch.reshape(len(BLM_2dsynch), 1) # so we can use nnls
    y = interpolate.interp1d(*aggregate_fill.blm_ir3())(aligned_secs)
    coef = nnls(BLM_2dsynch, y)[0]
    return coef*BLM_2dsynch

def halign(secs, losses, aggr_fill):
    """ Aligns losses to aggr_fill by returning a new 'secs' array """

    vloss_peak, iloss_peak = oml.imax(losses)
    vfill_peak, ifill_peak = oml.imax(aggr_fill.blm_ir3().y)
    delta = secs[iloss_peak] - aggr_fill.blm_ir3().x[ifill_peak]
    print("peaks\n\t2d-synch : {:.2f}\n\taggregate: {:.2f}\n\tdelta    : {:.2f}"
            .format(secs[iloss_peak], aggr_fill.blm_ir3().x[ifill_peak], delta))
    return secs - delta

def optimize(secs, losses, aggr_fill):
    """ Lets do a naive implementation:
            1. Remove a random particle
            2. Recalculate the lossmap. If one timestamp has lossmap < lhcdata, put back the particle
            3. If we have failed to remove a particle n number of times, stop. Else, continue with #1.
    """

    # We want aggr_fill to use the 'secs' timescale
    oml_y = interpolate.interp1d(*aggr_fill.blm_ir3())(secs)
    avg_loss = oml.moving_average(losses, BLM_INT)
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

        rstart = int(r - settings.BLM_INT/2)
        rend = int(r + settings.BLM_INT/2)
        temp_avg_loss[rstart:rend] -= 1.0/settings.BLM_INT

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


if __name__ == "__main__":
    ps = PhaseSpace(settings.STARTDIST_PATH)
    lossmap = lm.get_lossmap(settings.COLL_PATH)
    compare_to_LHC_aggregate(ps, lossmap)
