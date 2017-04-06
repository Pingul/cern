
import matplotlib.pyplot as plt
import numpy as np

from scipy import interpolate
from scipy.optimize import nnls 
from sklearn import linear_model
from bisect import bisect_left

import lossmap as lm
from phasespace import PhaseSpace
from settings import *

import sys, os
sys.path.append("/Users/swretbor/Workspace/collimation/proj/lhcstat")
import oml
import analyse_fill as af

sys.path.append("/Users/swretbor/Workspace/collimation/proj/common")
from logger import ModuleLogger, LogLevel

lg = ModuleLogger("lhccomp")

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

    min_turn = min(tm_lossmap.keys()) - settings.BLM_INT
    if prune:
        lg.log("pruning time scale")
        max_turn = int(16.5*11245.0)
        raise Exception("This has not been updates since all the changes -- probably invalid")
    else:
        lg.log("using full time scale")
        max_turn = max(tm_lossmap.keys()) + settings.BLM_INT

    lg.log("extract losses")
    turns = np.arange(min_turn, max_turn)
    secs = turns/11245.0
    losses = np.zeros(len(turns))
    for turn in tm_lossmap:
        losses[turn - min_turn] = len(tm_lossmap[turn])

    lg.log("emulate 2dsynch BLM")
    BLM_2dsynch = trailing_integration(losses, settings.BLM_INT)

    lg.log("align BLM")
    aggr_fill = af.aggregate_fill(beam, from_cache=True)
    asecs = halign(secs, BLM_2dsynch, aggr_fill)
    BLM_2dsynch = valign(aggr_fill, asecs, BLM_2dsynch)

    ## SEPARATE LOSSMAP
    lg.log("separate lossmap")
    lossmaps, action_values = lm.separate_lossmap(tm_lossmap, ps)
    action_values = np.array(action_values) - round(settings.H_SEPARATRIX)
    x = np.zeros(len(lossmaps)*len(turns), dtype=float).reshape(len(lossmaps), len(turns))
    for i, lossmap in enumerate(lossmaps):
        for turn in lossmap:
            x[i][turn - min_turn] = len(lossmap[turn])

        if integration:
            x[i] = trailing_integration(x[i], settings.BLM_INT)

    ## FIT
    lg.log("fit")
    xt = np.transpose(x)
    y = interpolate.interp1d(*aggr_fill.blm_ir3())(asecs)
    coef = nnls(xt, y)[0]
    coef_file = "lhc_fit_coefficients.txt"
    with open(coef_file, "w") as f:
        f.write("Coefficients:")
        for i, c in enumerate(coef):
            f.write("{:>5} = {:<20.10f} (H = {})\n".format("a{}".format(i), c, action_values[i]))
    lg.log("saved coefficients to '{}'".format(coef_file))
    x_fit = np.sum(coef*xt, axis=1)


    if not integration:
        x_fit = trailing_integration(x_fit, settings.BLM_INT)
        x_fit = x_fit.reshape(len(x_fit), 1)
        correction = nnls(x_fit, y)[0]
        x_fit *= correction

    lg.log("fitting completed")
    lg.log("plot")

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

    shift = 2*11245
    vloss_peak, iloss_peak = oml.imax(losses[shift:])
    iloss_peak += shift

    vfill_peak, ifill_peak = oml.imax(aggr_fill.blm_ir3().y)
    delta = secs[iloss_peak] - aggr_fill.blm_ir3().x[ifill_peak]
    lg.log("peaks\n\t2d-synch : {:.2f}\n\taggregate: {:.2f}\n\tdelta    : {:.2f}"
            .format(secs[iloss_peak], aggr_fill.blm_ir3().x[ifill_peak], delta))
    return secs - delta


if __name__ == "__main__":
    ps = PhaseSpace(settings.STARTDIST_PATH)
    lossmap = lm.get_lossmap(settings.COLL_PATH)
    compare_to_LHC_aggregate(ps, lossmap)
