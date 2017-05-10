
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

from scipy import interpolate 
from scipy.optimize import nnls, least_squares
from sklearn import linear_model
from bisect import bisect_left

import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../lhcstat/"))

from settings import *
import oml
import analyse_fill as af
import lossmap as lm

from logger import ModuleLogger, LogLevel
from phasespace import PhaseSpace

lg = ModuleLogger("lhccomp")

FIT_COEF_FILE = "fit_coef.txt"

def trailing_integration(sequence, N):
    """ For each entry, integrate the N first entries.
    """

    sums = np.convolve(sequence, np.ones((N,)))
    return sums[:len(sequence)]


class LHCComparison:
    """ Compare a fill with data from the toy model """

    def __init__(self, fill, ps, lossmap):
        self.fill = fill
        self.ps = ps
        self.lossmap = lossmap

        self.min_t = min(self.lossmap.keys()) - settings.BLM_INT
        self.max_t = max(self.lossmap.keys()) + settings.BLM_INT

        self.secs = np.arange(self.min_t, self.max_t)/11245.0
        self.losses = np.zeros(self.secs.size)
        for turn in self.lossmap:
            self.losses[turn - self.min_t] = len(self.lossmap[turn])

        self.opt_mask = np.ones(self.secs.shape, dtype=bool)

    def set_window(self, t_start=None, t_end=None):
        lg.log("pruning time scale to {}-{}".format(t_start, t_end), log_level=LogLevel.notify)
        if not t_start is None:
            self.opt_mask *= self.secs > t_start
        if not t_end is None:
            self.opt_mask *= self.secs < t_end

    def halign(self):
        lg.log("manual halign", log_level=LogLevel.warning)
        self.secs -= 0.45 # this seems to be a good choice with no betatron amplitude

    def BLM(self, normalised=True):
        blm = trailing_integration(self.losses, settings.BLM_INT)
        if (normalised):
            fmax = self.fill.blm_ir3().y.max()
            blm /= blm.max()/fmax
        return blm[self.opt_mask]

    def t(self):
        """ Timescale """
        return self.secs[self.opt_mask]

    def fit_action_values(self):
        lms, avs = lm.separate_lossmap(self.lossmap, self.ps, True)
        avs = np.array(avs)
        for i in range(avs.size):
            if avs[i] < 0: avs[i] += round(settings.H_SEPARATRIX)
            else: avs[i] -= round(settings.H_SEPARATRIX)

        lg.log("separate lossmaps")
        blms = np.zeros((len(lms), self.secs.size), dtype=float)
        for i, l in enumerate(lms):
            for turn in l:
                blms[i][turn - self.min_t] = len(l[turn])

            if settings.TRAILING_INTEGRATION:
                blms[i] = trailing_integration(blms[i], settings.BLM_INT)
        blms = blms[:, self.opt_mask]

        lg.log("fit")
        blms = np.transpose(blms)
        y = interpolate.interp1d(*self.fill.blm_ir3())(self.t())

        self.action_values = avs
        self.coef = nnls(blms, y)[0]
        self.blm_fit = np.sum(self.coef*blms, axis=1)

    def fit_results(self):
        return { 
                 "blm_fit" : self.blm_fit,
                 "c" : self.coef,
                 "action_values" : self.action_values
                }


def plot_comp(fill, blm=None, fit=None, block=True):
    fig, loss_ax = plt.subplots()
    loss_ax.plot(*fill.blm_ir3(), label="Aggr. fill (beam {})".format(fill.beam), color='r')
    loss_ax.set_yscale("log")
    loss_ax.set_xlim([-5, 40])
    loss_ax.set_ylim([0.5e-5, 1])
    loss_ax.set_ylabel("Losses (∆particles/1.3s)")
    loss_ax.set_xlabel("t (s)")
    loss_ax.axvspan(0.0, fill.crossover_point()['t'], facecolor='b', zorder=0, alpha=0.1)

    if not blm is None:
        loss_ax.plot(*blm, label="Toy model BLM")
    if not fit is None:
        loss_ax.plot(*fit, label="Fit", linestyle='--', zorder=6, color="forestgreen")

    loss_ax.legend(loc="upper right")

    plt.title("Compare 2d-synchrotron with aggregate fill")
    if block:
        plt.show()
    else:
        plt.draw()



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
    aggr_fill = af.aggregate_fill(beam, from_cache=True)

    min_turn = min(tm_lossmap.keys()) - settings.BLM_INT
    max_turn = max(tm_lossmap.keys()) + settings.BLM_INT
    # if prune:
        # lg.log("pruning time scale", log_level=LogLevel.notify)
        # max_turn = int(aggr_fill.crossover_point()['t']*11245.0)
        # tm_lossmap = {k : tm_lossmap[k] for k in tm_lossmap if k < max_turn}
    # else:
        # lg.log("using full time scale")
        # max_turn = max(tm_lossmap.keys()) + settings.BLM_INT

    lg.log("extract losses")
    turns = np.arange(min_turn, max_turn)
    secs = turns/11245.0
    losses = np.zeros(len(turns))
    for turn in tm_lossmap:
        losses[turn - min_turn] = len(tm_lossmap[turn])

    lg.log("emulate 2dsynch BLM")
    BLM_2dsynch = trailing_integration(losses, settings.BLM_INT)

    ## SEPARATE LOSSMAP
    lg.log("separate lossmap")
    lossmaps, action_values = lm.separate_lossmap(tm_lossmap, ps)
    action_values = np.array(action_values) - round(settings.H_SEPARATRIX)

    # lg.log("temporary pruning action values", log_level=LogLevel.warning)
    # to_delete = np.where(action_values > 0)
    # lossmaps = np.delete(lossmaps, to_delete)
    # action_values = np.delete(action_values, to_delete)
    
    x = np.zeros(len(lossmaps)*len(turns), dtype=float).reshape(len(lossmaps), len(turns))
    for i, lossmap in enumerate(lossmaps):
        for turn in lossmap:
            x[i][turn - min_turn] = len(lossmap[turn])

        if settings.TRAILING_INTEGRATION:
            x[i] = trailing_integration(x[i], settings.BLM_INT)

    lg.log("align BLM")
    asecs = halign(secs, BLM_2dsynch, aggr_fill)
    BLM_2dsynch = valign(aggr_fill, asecs, BLM_2dsynch)

    if settings.PRUNE_TIMESCALE:
        a, b = (14.5, 20)
        lg.log("pruning time scale", log_level=LogLevel.warning)
        lg.log("optimising between {}-{}".format(a, b), log_level=LogLevel.notify)
        mask = asecs < b
        mask *= asecs > a
        # mask *= asecs < 14.5
        asecs = asecs[mask]
        x = x[:, mask]
        
    ## FIT
    lg.log("fit")
    xt = np.transpose(x)
    y = interpolate.interp1d(*aggr_fill.blm_ir3())(asecs)

    method = "linear"
    lg.log(method)
    if method == "log":
        xt[xt == 0] = 1e-5
        xt = np.log10(xt)
        # y = np.log10(y)
        coef = np.linalg.lstsq(xt, y)[0]
        x_fit = np.sum(coef*xt, axis=1)
    elif method == "linear":
        coef = nnls(xt, y)[0]
        x_fit = np.sum(coef*xt, axis=1)
    else:
        raise Exception("not valid")

    with open(FIT_COEF_FILE, "w") as f:
        f.write("# {:>3} {:<20} {}\n".format("n", "Coefficients", "∆H"))
        for i, c in enumerate(coef):
            f.write(" {:>3} {:<20.10f} {}\n".format(i, c, action_values[i]))
    lg.log("saved coefficients to '{}'".format(FIT_COEF_FILE))


    if not settings.TRAILING_INTEGRATION:
        x_fit = trailing_integration(x_fit, settings.BLM_INT)
        x_fit = x_fit.reshape(len(x_fit), 1)
        correction = nnls(x_fit, y)[0]
        x_fit *= correction

    lg.log("fitting completed")
    lg.log("plot")

    # Plotting lossmap fit
    option_string = "(integ.)" if settings.TRAILING_INTEGRATION else ""
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
    plot_coefficients(action_values, coef/coef.max(), plot_type="bar")

def valign(aggregate_fill, asecs, BLM_2dsynch):
    """ Returns a vertically shifted BLM_2dsynch array that is fitted 
        w.r.t the aggregate fill. We assume that they have been horizontally aligned
    """
    BLM_2dsynch = BLM_2dsynch.reshape(len(BLM_2dsynch), 1) # so we can use nnls
    y = interpolate.interp1d(*aggregate_fill.blm_ir3())(asecs)
    coef = nnls(BLM_2dsynch, y)[0]
    return coef*BLM_2dsynch

def halign(secs, losses, aggr_fill):
    """ Aligns losses to aggr_fill by returning a new 'secs' array """


    shift = 2*11245
    vloss_peak, iloss_peak = oml.imax(losses[shift:])
    iloss_peak += shift

    vfill_peak, ifill_peak = oml.imax(aggr_fill.blm_ir3().y)
    delta = secs[iloss_peak] - aggr_fill.blm_ir3().x[ifill_peak]
    lg.log("peaks\n\ttoymodel : {:.2f}\n\taggregate: {:.2f}\n\tdelta    : {:.2f}"
            .format(secs[iloss_peak], aggr_fill.blm_ir3().x[ifill_peak], delta))

    # lg.log("halign disabled", log_level=LogLevel.notify)
    # return secs

    # seems to work okay when there is no horizontal displacement
    lg.log("manual halign", log_level=LogLevel.warning)
    return secs - 0.43

    return secs - delta

def dist_curve(H, coef, ctype = "linear"):
    p_init = [0, 0, 0]
    if ctype == "linear":
        model = lambda x, k, m: k*x + m
    elif ctype == "exp_pow":
        # model = lambda x, k, m: np.exp(k*np.sqrt(x/m))
        # model = lambda x, k, m: np.exp(k*np.power(x/m, 0.5))
        model = lambda x, k, m, c: np.exp(k*np.power(x, m) + c)
        p_init = [-1.0, 1.0, 0]
    elif ctype == "exp":
        model = lambda x, k, m: np.exp(k*x + m)
    elif ctype == "log":
        model = lambda x, k, m: k*np.log(x) + m
    elif ctype == "1/x":
        model = lambda x, k, m: k/(x + m)
    else:
        raise Exception("not valid model type")

    func = lambda k, h, c: (model(h, *k) - c)**2
    res = least_squares(func, p_init, loss='soft_l1', f_scale=0.0001, args=(H, coef), jac='2-point')
    # res = least_squares(func, p_init, loss='linear', f_scale=0.05, args=(H, coef), jac='3-point')
    #                                                        ^^^^ seems to determine effect of outliers (lower, less influence)
    # lg.log(res)
    lg.log(res.message)
    lg.log(res.optimality)
    lg.log("Fit: ", *zip(("k", "m", "c"), res.x))
    return lambda x : model(x, *res.x)

def plot_coefficients(H, coef, plot_type="bar", curve=None, info=None, block=True):
    title = "Optimised action value coefficients"
    if not info is None:
        title += ": {}".format(info)

    if plot_type == "scatter":
        fig, ax = plt.subplots()
        if curve:
            ax.plot(H, list(map(curve, H)), color="darkorange", label='fit')
        i = 0
        ax.scatter(H[i:], coef[i:], label='used')
        ax.scatter(H[:i], coef[:i], color='gray', label='unused')
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed')
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.legend(loc='upper right')

        ax.set_ylabel("Value")
        ax.set_xlabel("∆H")

        plt.title(title)
        plt.tight_layout()

    elif plot_type == "bar":
        fig, ax = plt.subplots()
        ax.bar(range(len(coef)), coef)
        ax.set_xticks(range(len(coef)))
        ax.set_xticklabels(["{:.3g}".format(h) for h in H], rotation='vertical')
        index = bisect_left(H, 0) - 0.5
        ax.axvspan(index - 0.025, index + 0.025, facecolor='r', zorder=0, alpha=0.2)
        ax.text(index, max(coef)/2.0, "Separatrix", fontsize=10, 
                ha='center', va='center', rotation='vertical', color='r')

        ax.set_ylabel("Value")
        ax.set_xlabel("∆H")

        plt.title(title)
        plt.tight_layout()
    elif plot_type == "scale":
        fig = plt.figure()
        ax = fig.add_subplot(221)
        ax.scatter(H, coef)
        ax.yaxis.grid(color='gray', linestyle='dashed')
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.set_ylabel("Value")
        ax.set_xlabel("∆H")

        ax2 = fig.add_subplot(222)
        ax2.scatter(H, coef)
        ax2.set_xscale("log")
        ax2.yaxis.grid(color='gray', linestyle='dashed')
        ax2.xaxis.grid(color='gray', linestyle='dashed')
        ax2.set_ylabel("Value")
        ax2.set_xlabel("∆H [log]")

        ax3 = fig.add_subplot(223)
        ax3.scatter(H, coef)
        ax3.set_yscale("log")
        ax3.yaxis.grid(color='gray', linestyle='dashed')
        ax3.xaxis.grid(color='gray', linestyle='dashed')
        ax3.set_ylabel("Value [log]")
        ax3.set_xlabel("∆H")

        ax4 = fig.add_subplot(224)
        ax4.scatter(H, coef)
        ax4.set_yscale("log")
        ax4.set_xscale("log")
        ax4.yaxis.grid(color='gray', linestyle='dashed')
        ax4.xaxis.grid(color='gray', linestyle='dashed')
        ax4.set_ylabel("Value [log]")
        ax4.set_xlabel("∆H [log]")
        fig.suptitle("Coefficient scale plot {}".format(info if not info is None else ""))

        if (curve):
            ax.plot(H, list(map(curve, H)), color="darkorange")
            ax2.plot(H, list(map(curve, H)), color="darkorange")
            ax3.plot(H, list(map(curve, H)), color="darkorange")
            ax4.plot(H, list(map(curve, H)), color="darkorange")
    else:
        raise Exception("plot_type not valid")

    if block:
        plt.show()
    else:
        plt.draw()

if __name__ == "__main__":
    action = sys.argv[1]
    
    fill = af.aggregate_fill(1, from_cache=True)
    ps = PhaseSpace(settings.STARTDIST_PATH)
    lossmap = lm.get_lossmap(settings.COLL_PATH)
    comp = LHCComparison(fill, ps, lossmap)

    if action == "fit":
        blm_s = (comp.t(), comp.BLM())
        # comp.set_window(5, 11)
        comp.halign()
        comp.fit_action_values()

        r = comp.fit_results()
        coef = r['c']
        avs = r['action_values']
        with open(FIT_COEF_FILE, "w") as f:
            f.write("# {:>3} {:<20} {}\n".format("n", "Coefficients", "∆H"))
            for i, c in enumerate(coef):
                f.write(" {:>3} {:<20.16f} {}\n".format(i, c, avs[i]))
        lg.log("saved coefficients to '{}'".format(FIT_COEF_FILE))

        blm_fit = (comp.t(), r['blm_fit'])
        plot_comp(fill, blm_s, blm_fit, block=False)
        plot_coefficients(avs, coef)
    elif action == "compare":
        plot_comp(fill, (comp.t(), comp.BLM()))
    elif action == "fit_curve":
        coef, H = np.loadtxt(FIT_COEF_FILE, skiprows=1, usecols=(1, 2), unpack=True)
        if np.any(H < 0):
            # assuming outside
            lg.log("below", log_level=LogLevel.notify)
            below = np.where(H < 0)
            hb = -H[below]
            cb = coef[below]
            
            cb /= cb.max()
            hb /= hb.max()

            print(cb.size)
            m = np.invert(cb < 1e-9)
            m *= np.invert((hb < 0.3)*(cb < 5e-3))
            # m = np.invert((hb < 0.05)*(cb < 1e-5))
            # m *= np.invert((hb > 0.2)*(hb < 0.4)*(cb < 2e-3))
            m *= np.invert((hb < 0.2)*(cb < 1e-1))
            m *= np.invert((hb < 0.3)*(cb < 2e-2))
            m *= np.invert((hb < 0.1)*(cb < 0.25))
            hb = hb[m]
            cb = cb[m]
            print(cb.size)

            f = dist_curve(hb, cb, "exp_pow")
            plot_coefficients(hb, cb, plot_type="scale", curve=f)
        else:
            coef /= coef.max()
    
            # outside
            # f = dist_curve(H, coef, "exp_pow")
            # H /= H.max()
            # plot_coefficients(H, coef, plot_type="scale", curve=f)
    
            # inside
            mask = coef > 0
            H = H[mask]
            coef = coef[mask]
            f = dist_curve(H, coef, "linear")
            plot_coefficients(H, coef, plot_type="scatter", curve=f)

    elif action == "test":
        fill = af.aggregate_fill(1, from_cache=True)
        ps = PhaseSpace(settings.STARTDIST_PATH)
        lossmap = lm.get_lossmap(settings.COLL_PATH)

        comp = LHCComparison(fill, ps, lossmap)
        blm_s = (comp.t(), comp.BLM())
        comp.set_window(10, 20)
        comp.fit_action_values()
        r = comp.fit_results()
        blm_fit = (comp.t(), r['blm_fit'])
        plot_comp(fill, blm_s, blm_fit, block=False)
        plot_coefficients(r['action_values'], r['c'])
    else:
        lg.log("unrecognised action", log_level=LogLevel.warning)
