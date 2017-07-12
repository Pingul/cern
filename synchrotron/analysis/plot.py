# coding=utf-8
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import itertools
import glob

from sys import argv

from phasespace import PhaseSpace
from lossmap import *
from settings import *

from logger import ModuleLogger, LogLevel
lg = ModuleLogger("plot")

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
                # lg.log(j, id)
                series[key]['x'][i*nbr_p + j] = i
                series[key]['y'][i*nbr_p + j] = ps.h[i*ps.nbr_p + id]
    
    ax.scatter(series['alive']['x'][:-1], series['alive']['y'][:-1], color='b', s=1, zorder=10)
    ax.scatter(series['lost']['x'][:-1], series['lost']['y'][:-1], color='r', s=4)


    # Should in theory do the same thing as the code above, but it's way slower
    # x = np.empty(ps.h.shape)
    # c = np.empty(ps.h.shape)
    # for i in range(ps.nbr_turns):
        # for j in range(ps.nbr_p):
            # x[i*ps.nbr_p + j] = i
            # c[i*ps.nbr_p + j] = 1 if j in pbin['alive'] else 100

    # lg.log("drawing {} points...".format(len(x)))
    # ax.scatter(x, ps.h, c=c, s=4)

    ax.set_xlabel("Saved turn")
    ax.set_ylabel("Hamiltonian")
    plt.show()

def plot_hamiltonian_dist_histogram(ps):
    fig, ax = plt.subplots()
    h_val = np.array(ps.h - settings.H_SEPARATRIX, dtype=int)
    hmax = h_val.max()
    hmin = h_val.min()
    nbr_bins = 100
    bins = np.arange(hmin, hmax, (hmax-hmin)/nbr_bins)
    ax.hist(h_val, bins=bins, edgecolor='white')
    ax.set_title("Action value H starting distribution")
    ax.set_xlabel("∆H")
    ax.set_ylabel("#")
    plt.show()

def plot_e_dist_histogram(ps):
    fig, ax = plt.subplots()
    nbr_bins = 150
    emin = ps.denergy.min()
    emax = ps.denergy.max()
    bins = np.arange(emin, emax, (emax-emin)/nbr_bins)
    ax.hist(ps.denergy, bins=bins, edgecolor='white')
    ax.set_title("Action value ∆E distribution")
    ax.set_xlabel("∆E")
    ax.set_ylabel("#")
    plt.show()

phasespace_dir = "/Users/swretbor/Workspace/work_afs/2dsynch/phasespace/phasespace"
def phasespace_frame(num, ps):
    filename = "{}/{}lines.dat".format(phasespace_dir, num)
    ps.plot_trajectory(filename, redraw=True)

def phasespace_evolution():
    frames = 300;

    ps = PhaseSpace(None)
    ps.create_plot()
    ps.format_axes()
    ps.plot_trajectory("{}/0lines.dat".format(phasespace_dir))

    ani = animation.FuncAnimation(ps.fig, phasespace_frame, int(frames), fargs=(ps,), interval=100, blit=False)
    mov_file = "{}/mov.mp4".format(phasespace_dir) 
    lg.log("saving movie to '{}'".format(mov_file))
    ani.save(mov_file)
    lg.log("finished saving to '{}'".format(mov_file))


SAVE_FILE = ''
if (len(argv) > 2):
    SAVE_FILE = argv[2]

if __name__ == "__main__":
    ACTION = argv[1]
    if ACTION == "animate":
        lg.log("animate trajectory")
        ps = PhaseSpace(settings.PARTICLE_PATH)
        ps.animate(get_lossmap(settings.COLL_PATH), save_to=SAVE_FILE)
    if ACTION == "animate-full":
        lg.log("animating trajectory and background")
        output = "calc/animate-full.mp4"
        lg.log("will save movie to '{}'".format(output))
        ps = PhaseSpace(settings.PARTICLE_PATH)
        ps.animate(get_lossmap(settings.COLL_PATH), output, animateBackground=True)
    elif ACTION == "lossmap":
        lg.log("plot lossmap")
        lossmap = get_lossmap(settings.COLL_PATH)
        plot_lossmap([lossmap], ["Toy model"], SAVE_FILE)
    elif ACTION == "lossmap2":
        lg.log("plot lossmap2")
        ps = PhaseSpace(settings.STARTDIST_PATH)
        lms = []
        labels = []
        for filepath in glob.glob(settings.DATA_DIR + "/*.ch"):
            lms.append(get_lossmap(filepath))
            labels.append(filepath.split('/')[-1].split('.')[0])
        plot_lossmap(lms, labels)

    elif ACTION == "separated-lossmap":
        lg.log("plot one series for each action value of the lossmap")
        ps = PhaseSpace(settings.STARTDIST_PATH)
        lossmap = get_lossmap(settings.COLL_PATH)
        lossmaps, labels = separate_lossmap(lossmap, ps, False)
        # plot_lossmap(lossmaps[0:-1:10], labels[0:-1:10])
        plot_lossmap(lossmaps, labels)
    elif ACTION == "energy":
        lg.log("plot energy oscillations")
        plot_energy_oscillations()
    elif ACTION == "startdist":
        lg.log("plot start distribution")
        losses_from_sec = float(argv[2]) if len(argv) > 2 else 0.0
        ps = PhaseSpace(settings.STARTDIST_PATH)
        lossmap = get_lossmap(settings.COLL_PATH)
        pbin = ps.categorize_particles(lossmap, losses_from_sec)
        ps.plot_particles(pbin)
    elif ACTION == "dist":
        lg.log("plot dist")
        if not len(argv) > 2: raise Exception("need to give a file path to plot")
        ps = PhaseSpace(argv[2])
        ps.plot_particles()
    elif ACTION == "phasespace":
        lg.log("plot phase space")
        ps = PhaseSpace(pfile=None)
        ps.create_plot()
        ps.plot_background_lines()
        ps.format_axes()
        if SAVE_FILE:
            plt.savefig(SAVE_FILE)
        else:
            plt.show()
    elif ACTION == "phasespace-mov":
        lg.log("animating phasespace evolution")
        phasespace_evolution()
    elif ACTION == "trajectory":
        lg.log("plot particle trajectory")
        ps = PhaseSpace(settings.PARTICLE_PATH)
        ps.create_plot()
        ps.format_axes()
        ps.plot_background_lines()
        ps.plot_trajectory(randomizeColors=True)
        plt.show()
    elif ACTION == "ham-evolution":
        lg.log("plot hamiltonian")
        ps = PhaseSpace(settings.PARTICLE_PATH)
        plot_hamiltonian_evolution(ps)
    elif ACTION == "ham-dist":
        lg.log("plot action distribution")
        input_file = settings.STARTDIST_PATH
        if len(argv) == 3: input_file = argv[2]
        ps = PhaseSpace(input_file)
        plot_hamiltonian_dist_histogram(ps)
    elif ACTION == "e-dist":
        lg.log("plot energy distribution")
        input_file = settings.STARTDIST_PATH
        if len(argv) == 3: input_file = argv[2]
        ps = PhaseSpace(input_file)
        plot_e_dist_histogram(ps)

    elif ACTION == "x":
        lg.log("plotting x distribution")
        input_file = settings.STARTDIST_PATH
        if len(argv) == 3: input_file = argv[2]
        ps = PhaseSpace(input_file)
        fig, ax = plt.subplots()
        lm = get_lossmap(settings.COLL_PATH)
        pbin = ps.categorize_particles(lm)
        
        plot_first = True
        if plot_first:
            if len(pbin['alive']) > 0: ax.scatter(ps.x[pbin['alive']], ps.px[pbin['alive']]);
            if len(pbin['lost']) > 0: ax.scatter(ps.x[pbin['lost']], ps.px[pbin['lost']], color='r');
        else:
            ax.scatter(ps.x, ps.px)

        ax.set_xlabel("x")
        ax.set_ylabel("x'")
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray')
        ax.xaxis.grid(color='gray')
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.1E}".format(x)))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.1E}".format(x)))
        plt.title("Motion for geometric coordinates at IR3 TCP")
        plt.show()
    elif ACTION == 'impacts':
        plot_first_impacts(settings.COLL_PATH)
    elif ACTION == "fort13":
        input_file = argv[2]
        ps = PhaseSpace.from_fort13(input_file)
        ps.plot_particles()
    elif ACTION == 'dist0':
        input_file = argv[2]
        ps = PhaseSpace.from_dist0(input_file)
        ps.plot_particles()
    elif ACTION == 'dump':
        input_file = argv[2]
        ps = PhaseSpace.from_six_dump(input_file)
        ps.plot_particles()
    else:
        lg.log("unrecognised action")
