# coding=utf-8
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import itertools

from sys import argv

from phasespace import PhaseSpace
from lossmap import *
from settings import *


def plot_energy_oscillations():
    nbr_turns = 500*11245
    ramp = read_ramp(settings.RAMP_PATH, nbr_turns)
    turns = np.array(range(nbr_turns))

    fig = plt.figure()
    e_ax = fig.add_subplot(311)
    de_ax = fig.add_subplot(312, sharex=e_ax)
    ph_ax = fig.add_subplot(313, sharex=e_ax)

    e_ax.plot(turns, ramp['e'], color='b')
    # e_ax.plot(turns, ramp['e_fitted'], color='red')
    e_ax.set_ylabel("E (GeV)")
    e_ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1000.0)))
    e_ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e9)))

    de_ax.plot(turns, ramp['de'], color='b')
    # de_ax.plot(turns, ramp['de_fitted'], color='red')
    de_ax.set_ylabel("∆E (MeV/turn)")
    de_ax.set_xlabel("Time (kturns)")
    de_ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e6)))

    k = 2.9491187074838457087e-07
    phase = np.arcsin(ramp['de']/((k*turns + 6)*1e6))
    ph_ax.plot(turns, phase, color='gray')
    ph_ax.set_ylabel("φ_s (rad)")

    fig.suptitle("LHC ramp")

    plt.show()

def plot_hamiltonian_evolution(ps): # input phase space containing ps.h
    """ ps : PhaseSpace
            Should be time dependent (e.g. 'ps = PhaseSpace(settings.PARTICLE_PATH)')
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    lossmap = get_lossmap(settings.COLL_PATH)
    pbin = ps.categorize_particles(lossmap)
    series = {}
    for key in pbin: # 'alive'/'lost'
        nbr_p = len(pbin[key])
        size = ps.nbr_turns*nbr_p
        series[key] = {'x' : np.empty( (size,) ), 'y' : np.empty( (size,) )}
        for i in range(ps.nbr_turns):
            for j, id in enumerate(pbin[key]):
                series[key]['x'][i*nbr_p + j] = i
                series[key]['y'][i*nbr_p + j] = ps.h[i*ps.nbr_p + id]

    ax.scatter(series['alive']['x'], series['alive']['y'], color='b', s=1, zorder=10)
    ax.scatter(series['lost']['x'], series['lost']['y'], color='r', s=4)


    # Should in theory do the same thing as the code above, but it's way slower
    # x = np.empty(ps.h.shape)
    # c = np.empty(ps.h.shape)
    # for i in range(ps.nbr_turns):
        # for j in range(ps.nbr_p):
            # x[i*ps.nbr_p + j] = i
            # c[i*ps.nbr_p + j] = 1 if j in pbin['alive'] else 100

    # print("drawing {} points...".format(len(x)))
    # ax.scatter(x, ps.h, c=c, s=4)

    ax.set_xlabel("Saved turn")
    ax.set_ylabel("Hamiltonian")
    plt.show()

def plot_hamiltonian_dist_histogram(ps):
    fig, ax = plt.subplots()
    h_val = np.array(ps.h - settings.H_SEPARATRIX, dtype=int)
    hmax = h_val.max()
    hmin = h_val.min()
    nbr_bins = 40
    bins = np.arange(hmin, hmax, (hmax-hmin)/nbr_bins)
    ax.hist(h_val, bins=bins, edgecolor='white')
    ax.set_title("Action value H starting distribution")
    ax.set_xlabel("∆H")
    ax.set_ylabel("#")
    plt.show()


def plot_hamiltonian_lost_dist(ps, lm, time=16.5):
    """ ps : PhaseSpace
            Should be the start distribution. 
        lm : lossmap
    """
    pbin = ps.categorize_particles(lm, time)
    fig, ax = plt.subplots()
    data = {}
    for key in sorted(pbin): # sorted to make sure we get the same colors all the time
        data[key] = {}
        for pid in pbin[key]:
            h = round(ps.h[pid] - H_SEPARATRIX)
            if h in data[key]: data[key][h] += 1
            else: data[key][h] = 1
        
        if key == 'lost':
            label = "lost after {} s".format(time)
        elif key == 'discarded':
            label = "lost before {} s".format(time)
        else:
            label = key
        ax.scatter(list(data[key].keys()), list(data[key].values()), label=label)
    ax.legend(loc='upper left')
    ax.set_ylabel("Particles")
    ax.set_xlabel("∆H")
    plt.title("Hamiltonian lost distribution")
    plt.show()


def plot_lost():
    # lossmap = get_lossmap(settings.COLL_PATH, ["id"])
    df = pd.read_csv("calc/lost.dat", names=["id", "turn_lost", "coll_hit", "delta"])
    fig, ax = plt.subplots()
    ax.scatter(df["turn_lost"], df["delta"], s=4, color="b", label="travel time to collimator")
    ax.scatter(df["turn_lost"], df["coll_hit"], s=4, color="r", label="collimator hit")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.2f}".format(x/11245.0)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.2f}".format(x/11245.0)))
    ax.set_xlabel("t (s) when found outside of bucket")
    ax.set_ylabel("t (s)")
    ax.legend(loc="upper left")
    plt.show()


def phasespace_frame(num, ps):
    filename = "phasespace/{}lines.dat".format(num)
    ps.plot_trajectory(filename, redraw=True)

def phasespace_evolution():
    frames = 300;

    ps = PhaseSpace(None)
    ps.create_plot()
    ps.format_axes()
    ps.plot_trajectory("phasespace/0lines.dat")

    ani = animation.FuncAnimation(ps.fig, phasespace_frame, int(frames), fargs=(ps,), interval=100, blit=False)
    mov_file = "phasespace/mov.mp4" 
    print("saving movie to '{}'".format(mov_file))
    ani.save(mov_file)
    print("finished saving to '{}'".format(mov_file))


SAVE_FILE = ''
if (len(argv) > 2):
    SAVE_FILE = argv[2]

if __name__ == "__main__":
    ACTION = argv[1]
    if ACTION == "animate":
        print("animate trajectory")
        ps = PhaseSpace(settings.PARTICLE_PATH)
        ps.animate(get_lossmap(settings.COLL_PATH), save_to=SAVE_FILE)
    if ACTION == "animate-full":
        print("animating trajectory and background")
        output = "calc/animate-full.mp4"
        print("will save movie to '{}'".format(output))
        ps = PhaseSpace(settings.PARTICLE_PATH)
        ps.animate(get_lossmap(settings.COLL_PATH), output, animateBackground=True)
    elif ACTION == "lossmap":
        print("plot lossmap")
        lossmap = get_lossmap(settings.COLL_PATH)
        plot_lossmap([lossmap], ["Toy model"], SAVE_FILE)
    elif ACTION == "separated-lossmap":
        print("plot one series for each action value of the lossmap")
        ps = PhaseSpace(settings.STARTDIST_PATH)
        lossmap = get_lossmap(settings.COLL_PATH)
        lossmaps, labels = separate_lossmap(lossmap, ps)
        plot_lossmap(lossmaps, labels)
    elif ACTION == "energy":
        print("plot energy oscillations")
        plot_energy_oscillations()
    elif ACTION == "startdist":
        print("plot start distribution")
        losses_from_sec = float(argv[2]) if len(argv) > 2 else 0.0
        ps = PhaseSpace(settings.STARTDIST_PATH)
        lossmap = get_lossmap(settings.COLL_PATH)
        pbin = ps.categorize_particles(lossmap, losses_from_sec)
        ps.plot_particles(pbin)
    elif ACTION == "dist":
        print("plot dist")
        if not len(argv) > 2: raise Exception("need to give a file path to plot")
        ps = PhaseSpace(argv[2])
        ps.plot_particles()
    elif ACTION == "phasespace":
        print("plot phase space")
        ps = PhaseSpace(pfile=None)
        ps.create_plot()
        ps.plot_background_lines()
        ps.format_axes()
        if SAVE_FILE:
            plt.savefig(SAVE_FILE)
        else:
            plt.show()
    elif ACTION == "phasespace-mov":
        print("animating phasespace evolution")
        phasespace_evolution()
    elif ACTION == "trajectory":
        print("plot particle trajectory")
        ps = PhaseSpace(settings.PARTICLE_PATH)
        ps.create_plot()
        ps.format_axes()
        ps.plot_background_lines()
        ps.plot_trajectory(randomizeColors=True)
        plt.show()
    elif ACTION == "ham-evolution":
        print("plot hamiltonian")
        ps = PhaseSpace(settings.PARTICLE_PATH)
        plot_hamiltonian_evolution(ps)
    elif ACTION == "ham-lostdist":
        print("ploting hamiltonian lost distribution")
        ps = PhaseSpace(settings.STARTDIST_PATH)
        lm = get_lossmap(settings.COLL_PATH)
        plot_hamiltonian_lost_dist(ps, lm, 16.5)
    elif ACTION == "lost":
        print("lost plot")
        plot_lost()
    elif ACTION == "ham-dist":
        print("plot action distribution")
        input_file = settings.STARTDIST_PATH
        if len(argv) == 3: input_file = argv[2]
        ps = PhaseSpace(input_file)
        plot_hamiltonian_dist_histogram(ps)
    else:
        print("unrecognised action")