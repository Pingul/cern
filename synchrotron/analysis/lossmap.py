import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import numpy as np

from scipy import stats
from settings import settings
from ramp import read_ramp

import sys
sys.path.append("../../common")
from logger import ModuleLogger, LogLevel

lg = ModuleLogger("lossmap")

def trailing_integration(sequence, N):
    """ For each entry, integrate the N first entries.
    """

    sums = np.convolve(sequence, np.ones((N,)))
    return sums[:len(sequence)]

class CHitMap:
    """ Collimator hit map """
    
    def __init__(self, collfile="", pid_offset=0):
        if collfile != "":
            with open(collfile, 'r') as f:
                self.nbr_p = int(f.readline().strip())
                self.ids, self.turns, self.phase, self.denergy, self.x = np.loadtxt(f.readlines(), delimiter=',', skiprows=1, unpack=True)
                self.ids = self.ids.astype(int) + pid_offset
                self.turns = self.turns.astype(int)

    def old_lossmap(self):
        """ compatible with what 'get_lossmap' returned previously """ 
        hits = {}
        for i, pid in enumerate(self.ids):
            try: hits[self.turns[i]].append(pid)
            except: hits[self.turns[i]] = [pid]
        return hits

    def span(self):
        """ Returns the extent in time [turn] for first/last loss """
        return (self.turns.min(), self.turns.max())

    def losses(self, integrated=False):
        """ Returns 1d array with number of hits/turn """ 
        l = np.zeros(self.turns.max() + 1, dtype=int)
        for t in self.turns:
            l[t] = l[t] + 1

        if integrated:
            l = np.concatenate((l, np.zeros(settings.BLM_INT)))
            l = trailing_integration(l, settings.BLM_INT)
        return l

    def split(self, ps, separate_above_bucket=False):
        """ Separate the CHitMap into one hitamp per discrete action values, 
            given the initial PhaseSpace distribution 'ps'.

            Set separate_above_bucket=True if a distinction between particles ±∆E should be made.
        """
        h = ps.h.astype(int)
        uh = np.unique(h)
        ch = np.empty(uh.size, dtype=object)
        tot_p = 0
        for i, v in enumerate(uh):
            ch[i] = CHitMap()
            m = v == h
            ch[i].ids = self.ids[m]
            ch[i].turns = self.turns[m]
            ch[i].phase = self.phase[m]
            ch[i].denergy = self.denergy[m]
            ch[i].x = self.x[m]
            ch[i].nbr_p = np.sum(m)
            tot_p += ch[i].nbr_p

        # Adding the check here as a precaution
        if not tot_p == self.nbr_p:
            raise Exception("Number of particles was not conserved {} != {}".format(tot_p, self.nbr_p))
        return ch

def plot(hitmaps, labels=[], save_to=''):
    colors = 'r'
    if len(hitmaps) == 2: 
        colors = ('r', 'g')
    else: 
        colors = plt.cm.Set3(np.linspace(0, 1, len(hitmaps)))

    if len(labels) == 0:
        labels = [str(v) for v in range(len(hitmaps))]
    
    max_t = 0
    fig, ax = plt.subplots()
    for i, hm in enumerate(hitmaps):
        max_t = max(max_t, hm.span()[1])
        ax.plot(hm.losses(), label="{} (coll. hit)".format(labels[i]), alpha=0.5, c=colors[i])
        ax.plot(hm.losses(integrated=True), label="{} (integ.)".format(labels[i]), zorder=4, c=colors[i])
    ax.set_ylabel("Losses")
    ax.set_xlabel("t (s)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.2f}".format(x/11245.0)))
    ax.spines['right'].set_position(('axes', 1.15))
    ax.set_yscale('log')
    ax.legend(loc='upper right')

    # # RAMP
    ramp = read_ramp(settings.RAMP_PATH, max_t)
    e_ax = ax.twinx()
    e_ax.plot(ramp, color='gray', zorder=0, label='LHC energy ramp')
    e_ax.set_axis_off()

    title = "Toy model off-momentum losses"
    if save_to: title += " ('{}')".format(save_to)
    plt.title(title)

    if (save_to): 
        lg.log("saving plot to {}".format(save_to))
        plt.savefig(save_to) 
    plt.show()


def hits_from_collfile(collfile, pid_offset=0, mute=False):
    hits = []
    with open(collfile, 'r') as f:
        if not mute: lg.log("reading collimation file '{}'".format(collfile))
        for i, line in enumerate(f):
            if i < 2: continue
            line_c = line.rstrip().split(',')
            pid, turn = map(int, line_c[0:2])
            hits.append((turn, pid+pid_offset))
    return hits

def lossmap_from_hits(hits):
    hits.sort(key=lambda x: x[0])
    coll_hits = {} 
    for hit in hits:
        turn, pid = hit
        try: coll_hits[turn].append(pid)
        except: coll_hits[turn] = [pid]

    # converting to numpy arrays
    # for turn in coll_hits:
        # coll_hits[turn] = np.array(coll_hits[turn])

    return coll_hits


def get_lossmap(collfile):
    """ Returns data as
        lossmap = { turn : [particles_lost] }
    """

    lg.log("creating lossmap")
    hits = hits_from_collfile(collfile)
    return lossmap_from_hits(hits)

def separate_lossmap(lossmap, phasespace, separate_above_bucket=False):
    """ Separate the lossmap on the action values. The function will bin particles into 'bins' size 
        lossmaps.

        returns ([lossmaps], [action values])
    """
    if separate_above_bucket:
        lg.log("'separate_above_bucket' : action values with ∆E < 0 will have a - sign in front", log_level=LogLevel.notify)
    div_lossmaps = {}
    for turn in lossmap:
        for pid in lossmap[turn]:
            h = int(phasespace.h[pid])
            if separate_above_bucket and phasespace.denergy[pid] < 0:
                h = -h

            if not h in div_lossmaps: 
                div_lossmaps[h] = {}
            lossm = div_lossmaps[h]
            if not turn in lossm:
                lossm[turn] = []
            lossm[turn].append(pid)
    # lossmaps = [div_lossmaps[action] for action in div_lossmaps]
    a_values = sorted(div_lossmaps.keys())
    lossmaps = [div_lossmaps[l] for l in a_values]
    return (lossmaps, a_values)

def plot_lossmap(lossmaps, labels=[], save_to=''):
    """
        lossmaps = [lossmap1, lossmap2, ...] 
    """
    if len(lossmaps) == 0:
        raise Exception("no losses found")

    max_turn = 0
    for lm in lossmaps:
        max_turn = max(max(max_turn, max(lm.keys())), 20*11245) # at least 20 s
    turns = np.array(range(max_turn + 100))

    # Plotting
    fig, loss_ax = plt.subplots()

    # LOSSES
    if len(lossmaps) == 1:
        color_list = 'r'
    elif (len(lossmaps) == 2):
        color_list = ['r', 'g']
    else:
        color_list = plt.cm.Set3(np.linspace(0, 1, len(lossmaps)))

    for i, lm in enumerate(lossmaps):
        # lm = lossmaps[action]
        losses = np.array([len(lm[turn]) if turn in lm else 0 for turn in turns])
        avg_loss = trailing_integration(losses, int(1.3*11245))

        loss_ax.plot(turns, avg_loss, color=color_list[i], label=("{} (integ.)".format(labels[i]) if len(labels) > i else ""), zorder=4)
        loss_ax.plot(turns, losses, color=color_list[i], label=("{} (coll. hit)".format(labels[i] if len(labels) > i else "")), zorder=3, alpha=0.5)

        # loss_ax.plot(turns, losses, color='red', label=("{} (coll. hit)".format(labels[i] if len(labels) > i else "")), zorder=3)
        # loss_ax.plot(turns, avg_loss, color="blue", label=("{} (BLM)".format(labels[i]) if len(labels) > i else ""), zorder=2, linestyle='--', alpha=0.9)
    loss_ax.set_ylabel("Losses (∆particles)")

    # loss_ax.set_ylabel("Losses (∆particles/turn)")
    loss_ax.set_xlabel("t (s)")
    loss_ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.2f}".format(x/11245.0)))
    loss_ax.spines['right'].set_position(('axes', 1.15))
    loss_ax.set_yscale('log')
    if len(labels) > 0:
        loss_ax.legend(loc='upper right')

    # # RAMP
    ramp = np.array(read_ramp(settings.RAMP_PATH, len(turns)))
    e_ax = loss_ax.twinx()
    e_ax.plot(turns, ramp, color='gray', zorder=0, label='LHC energy ramp')
    e_ax.set_axis_off()

    title = "Toy model off-momentum losses"
    if save_to: title += " ('{}')".format(save_to)
    plt.title(title)

    if (save_to): 
        lg.log("saving plot to {}".format(save_to))
        plt.savefig(save_to) 
    plt.show()

def plot_first_impacts(collfile):
    x = []
    with open(collfile, 'r') as f:
        f.readline()
        f.readline()
        for l in f.readlines():
            id, turn, ph, de, _x = map(float, l.strip().split(','))
            if turn < 1000: continue
            x.append(_x)
    x = np.array(x)

    fig, ax = plt.subplots()
    mx = x.max()
    mn = x.min()
    nbr_bins = 50
    bins = np.arange(mn, mx, (mx - mn)/nbr_bins)
    ax.hist(x, bins=bins, edgecolor='white')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.3f}".format(x*1e3)))
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("#")
    plt.title("First impacts")
    plt.show()

# def read_ramp(file, nbr_turns):
    # e = np.empty(nbr_turns)
    # with open(file, 'r') as f:
        # for i, line in enumerate(f.readlines()):
            # if i >= nbr_turns: break
            # e[i] = float(line.rstrip().split()[1])*1e6

    # turns = range(nbr_turns)
    # #de = np.gradient(e)
    # de = np.diff(e) # size n - 1
    # de = np.append(de, de[-1])

    # slope, intercept, r_value, p_value, std_err = stats.linregress(turns, de)
    # de_fitted = [slope*turn + intercept for turn in turns]
    # e_fitted = []
    # s = 0
    # for v in de_fitted:
        # s += v
        # e_fitted.append(s + 450e9)

    # return {'e' : e, 'e_fitted' : e_fitted, 'de' : de, 'de_fitted' : de_fitted}


