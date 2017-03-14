import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import numpy as np

from scipy import stats
from settings import *

def get_lossmap(collfile):
    """ Returns data as
        lossmap = { turn : [particles_lost] }
    """

    print("creating lossmap")
    hits = []
    with open(collfile, 'r') as f:
        print("reading collimation file '{}'".format(collfile))
        for i, line in enumerate(f):
            if i < 2: 
                continue
            line_c = line.rstrip().split(',')
            pid, turn = map(int, line_c[0:2])
            hits.append([turn, pid])

    hits.sort(key=lambda x: x[0])
    coll_hits = {} 
    for hit in hits:
        turn, pid = hit
        try: coll_hits[turn].append(pid)
        except: coll_hits[turn] = [pid]
    return coll_hits

def separate_lossmap(lossmap, phasespace):
    """ Separate the lossmap on the action values. The function will bin particles into 'bins' size 
        lossmaps.

        returns ([lossmaps], [action values])
    """
    div_lossmaps = {}
    for turn in lossmap:
        for pid in lossmap[turn]:
            h = int(phasespace.h[pid]) 
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
        max_turn = max(max_turn, max(lm.keys()))
    turns = np.array(range(max_turn + 100))

    # Plotting
    fig, loss_ax = plt.subplots()

    # LOSSES
    if len(lossmaps) == 1:
        color_list = 'r'
    else:
        color_list = plt.cm.Set3(np.linspace(0, 1, len(lossmaps)))

    for i, lm in enumerate(lossmaps):
        # lm = lossmaps[action]
        losses = np.array([len(lm[turn]) if turn in lm else 0 for turn in turns])
        # avg_loss = moving_average(losses, int(1.3*11245))
        avg_loss = losses
        loss_ax.plot(turns, avg_loss, color=color_list[i], label=(labels[i] if len(labels) > i else ""), zorder=3, alpha=0.9)
    # loss_ax.set_ylabel("Losses (∆particles/1.3s)")
    loss_ax.set_ylabel("Losses (∆particles/turn)")
    loss_ax.set_xlabel("t (s)")
    loss_ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.2f}".format(x/11245.0)))
    loss_ax.spines['right'].set_position(('axes', 1.15))
    loss_ax.set_yscale('log')
    if len(labels) > 0:
        loss_ax.legend(loc='upper right')

    # # RAMP
    ramp = np.array(read_ramp(RAMP_FILE, len(turns))['e'])
    ramp = np.gradient(ramp)
    e_ax = loss_ax.twinx()
    e_ax.plot(turns, ramp, color='gray', zorder=0, label='LHC energy ramp')
    e_ax.set_axis_off()

    fig.suptitle("Toy model off-momentum lossmap")
    if (save_to): 
        print("saving plot to {}".format(save_to))
        plt.savefig(save_to) 
    plt.show()

def read_ramp(file, nbr_turns):
    e = np.empty(nbr_turns)
    with open(file, 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i >= nbr_turns: break
            e[i] = float(line.rstrip().split()[1])*1e6

    turns = range(nbr_turns)
    #de = np.gradient(e)
    de = np.diff(e) # size n - 1
    de = np.append(de, de[-1])

    slope, intercept, r_value, p_value, std_err = stats.linregress(turns, de)
    de_fitted = [slope*turn + intercept for turn in turns]
    e_fitted = []
    s = 0
    for v in de_fitted:
        s += v
        e_fitted.append(s + 450e9)

    return {'e' : e, 'e_fitted' : e_fitted, 'de' : de, 'de_fitted' : de_fitted}


