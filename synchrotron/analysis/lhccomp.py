
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

    def __init__(self, fill, ps, hitmap):
        self.fill = fill
        self.ps = ps
        self.hitmap = hitmap

        self.min_t, self.max_t = self.hitmap.trange(integrated=True)
        self.secs = np.arange(self.min_t, self.max_t)/11245.0
        self.opt_mask = np.ones(self.secs.shape, dtype=bool)

    def set_window(self, t_start=None, t_end=None):
        lg.log("pruning time scale to {}-{}".format(t_start, t_end), log_level=LogLevel.notify)
        if not t_start is None:
            self.opt_mask *= self.secs > t_start
        if not t_end is None:
            self.opt_mask *= self.secs < t_end

    def halign(self, secs):
        lg.log("manual halign", log_level=LogLevel.warning)
        self.secs -= secs 

    def BLM(self, normalised=True):
        blm = self.hitmap.losses(integrated=True)
        if (normalised):
            fmax = self.fill.blm_ir3().y.max()
            blm /= blm.max()/fmax
        return blm[self.opt_mask]

    def t(self):
        """ Timescale """
        return self.secs[self.opt_mask]

    def fit_action_values(self):
        lg.log("separate hit maps")
        chs, avs = self.hitmap.split(self.ps, separate_above_bucket=False)
        for i in range(avs.size):
            if avs[i] < 0: avs[i] += round(settings.H_SEPARATRIX)
            else: avs[i] -= round(settings.H_SEPARATRIX)

        lg.log("integrating")
        blms = np.zeros((len(chs), self.secs.size), dtype=float)
        for i, hm in enumerate(chs):
            m = hm.trange(integrated=True)[1]
            blms[i][:m] = hm.losses(integrated=True)
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
        loss_ax.plot(*blm, label="Simulated BLM")
    if not fit is None:
        loss_ax.plot(*fit, label="Fit", linestyle='--', zorder=6, color="forestgreen")

    loss_ax.legend(loc="upper right")

    plt.title("Compare simulation with aggregate fill")
    if block:
        plt.show()
    else:
        plt.draw()


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
    hitmap = lm.CHitMap(settings.COLL_PATH)
    comp = LHCComparison(fill, ps, hitmap)

    if action == "fit":
        blm_s = (comp.t(), comp.BLM())
        comp.set_window(14, 20)
        # comp.halign()
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
        # comp.halign()
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
    else:
        lg.log("unrecognised action", log_level=LogLevel.warning)
