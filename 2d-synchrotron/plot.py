# coding=utf-8
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cmx
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle
import math
import numpy as np
import pandas as pd
import itertools
from scipy import stats

from sys import argv

## Some global definitions
PARTICLE_FILE = "calc/particles.dat"
LINE_FILE = "calc/lines.dat"
COLL_FILE = "calc/coll.dat"
STARTDIST_FILE = "calc/startdist.dat"
ENDDIST_FILE = "calc/enddist.dat"
RAMP_FILE = "resources/LHC_ramp.dat"

PLOT_FRAME = {
    'x' : [-0.1*math.pi, 2.1*math.pi],
    'y' : [-0.4e9, 0.4e9]
    # 'x' : [-2*math.pi, 4*math.pi],
    # 'y' : [-2e9, 2e9]
}

SAVE_FILE = ''
if (len(argv) > 2):
    SAVE_FILE = argv[2]

def moving_average(sequence, N):
    """ Moving average given the sequence. Returns an list equal in length
    to the one given """

    average = np.convolve(sequence, np.ones((N,))/N, mode='same')
    return average

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

class PhaseSpace:
    """ Only works for data output by the 2d synchrotron """

    def __init__(self, pfile, rfile=None):
        self.nbr_p = 0
        self.nbr_turns = 0
        self.denergy = np.array([])
        self.phase = np.array([])
        # self.ref_energy = []
        self.ref_energy = None
        self.background_lines = []


        if pfile:
            print("reading in particles from '{}'".format(pfile))
            with open(pfile, 'r') as file:
                # This does not work... why?
                # self.nbr_p, self.nbr_turns = map(int, file.readline().rstrip('\n').split(','))
                # for i in range(self.nbr_turns):
                    # line = file.readline().strip()
                    # if not line: break # EOF
                    # e_ref = float(line.strip())
                    # self.ref_energy.append(e_ref)
                    # for j in range(self.nbr_p):
                        # line = file.readline().strip()
                        # denergy, phase = map(float, line.split(','))
                        # self.denergy.append(denergy)
                        # self.phase.append(phase)

                for i, line in enumerate(file.readlines()):
                    if i == 0:
                        self.nbr_p, self.nbr_turns = map(int, line.rstrip('\n').split(','))
                        self.denergy = np.empty(self.nbr_p*self.nbr_turns)
                        self.phase = np.empty(self.nbr_p*self.nbr_turns)
                        self.h = np.empty(self.nbr_p*self.nbr_turns)
                        continue
                    denergy, phase, h = map(float, line.rstrip('\n').split(','))
                    self.denergy[i - 1] = denergy
                    self.phase[i - 1] = phase
                    self.h[i - 1] = h

        if rfile:
            raise Exception("removed support for ramp file")
            # print("reading in energy from '{}'".format(rfile))
            # with open(rfile, 'r') as file:
                # for i, line in enumerate(file.readlines()):
                    # if i >= self.nbr_turns:
                        # break
                    # ref_e = float(line.split()[-1])
                    # self.ref_energy.append(ref_e)

    def create_plot(self):
        self.fig, self.ax = plt.subplots()

    def plot_trajectory(self, filePath=None, randomizeColors=False, redraw=False):
        print("plot trajectory '{}'".format(filePath))
        nbr_series = 0
        series = []
        if filePath is not None:
            print("reading file '{}'".format(filePath))
            with open(filePath, 'r') as f:
                for i, line in enumerate(f.readlines()):
                    if i == 0:
                        nbr_series = int(line.rstrip().split(',')[0])
                        series = [{"denergy" : [], "phase" : [], "h" : 0.0} for l in range(nbr_series)]
                        continue
                    denergy, phase, h = map(float, line.rstrip("\n").split(","))
                    series[i % nbr_series]["phase"].append(phase)
                    series[i % nbr_series]["denergy"].append(denergy)
                    series[i % nbr_series]['h'] = h
        else:
            print("using local data")
            nbr_series = self.nbr_p
            series = [{"denergy" : [], "phase" : [], "h" : 0.0} for l in range(nbr_series)]
            for i, de in enumerate(self.denergy):
                series[i % nbr_series]["phase"].append(self.phase[i])
                series[i % nbr_series]["denergy"].append(de)
                series[i % nbr_series]['h'] = self.h[i]

        if redraw and not len(self.background_lines) == len(series):
            raise Exception("new frame does not match series in old")

        cMap = plt.get_cmap('plasma_r')
        cNorm = colors.Normalize(vmin=-5e4, vmax=4e5)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cMap)
        print("plotting series")
        for i, trail in enumerate(series):
            colorVal = scalarMap.to_rgba(trail['h'])
            if randomizeColors:
                colorVal = np.random.rand(3,)
            if redraw:
                self.background_lines[i].set_xdata(trail["phase"])
                self.background_lines[i].set_ydata(trail["denergy"])
                # self.ax.set_xdata(trail["phase"])
                # self.ax.set_ydata(trail["denergy"])
            else:
                line, = self.ax.plot(trail["phase"], trail["denergy"], '-', color=colorVal, zorder=1)
                self.background_lines.append(line)

    def plot_background_lines(self):
        self.plot_trajectory(LINE_FILE)

    def plot_collimators(self):
        print("plot collimators")
        try:
            with open(COLL_FILE, 'r') as f:
                f.readline()
                second_line = f.readline()
                top_coll, bot_coll = map(float, second_line.rstrip().split(','))
                self.collimator = {'top' : top_coll, 'bot' : bot_coll}
        except:
            print("could not read '{}', will not plot".format(COLL_FILE))
        else:
            self.coll_hits = get_lossmap(COLL_FILE)
            self.ax.axhspan(PLOT_FRAME['y'][0], self.collimator['bot'], facecolor='red', alpha=0.1)
            self.ax.axhspan(self.collimator['top'], PLOT_FRAME['y'][1], facecolor='red', alpha=0.1)

    def format_axes(self):
        self.ax.set_xlim(PLOT_FRAME['x'])
        self.ax.set_ylim(PLOT_FRAME['y'])
        self.ax.set_xlabel("Phase (radians)")
        self.ax.set_ylabel("∆E (GeV)")
        self.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e9)))
        self.ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.1f}π".format(x/math.pi)))

    def categorize_particles(self, lossmap, sec=0.0):
        """ Seperates particle id's into two lists, 'lost' and 'alive' after time 'sec'.
            Particles are 'discarded' if they got lost before 'sec'. """
        pbin = {'lost' : [], 'discarded' : []}
        for turn in lossmap:
            if turn >= 11245.0*sec:
                pbin['lost'] += [loss for loss in lossmap[turn] if loss < self.nbr_p]
            else:
                pbin['discarded'] += [loss for loss in lossmap[turn] if loss < self.nbr_p]
        pbin['lost'].sort()
        pbin['discarded'].sort()
        pbin['alive'] = [i for i in range(self.nbr_p) if i not in pbin['lost'] and i not in pbin['discarded']]
        return pbin

    def plot_particles(self, pbin=None):
        self.create_plot()
        self.plot_background_lines()
        self.plot_collimators()
        self.format_axes()
        
        if pbin:
            print("lost:", len(pbin['lost']), 'alive:', len(pbin['alive']))
            live_z = 10
            lost_z = 9 if len(pbin['lost']) > len(pbin['alive']) else 11
            self.ax.scatter([self.phase[i] for i in pbin['alive']], [self.denergy[i] for i in pbin['alive']], color='b', s=2, zorder=live_z)
            self.ax.scatter([self.phase[i] for i in pbin['lost']], [self.denergy[i] for i in pbin['lost']], color='r', s=2, zorder=lost_z)
            self.ax.scatter([self.phase[i] for i in pbin['discarded']], [self.denergy[i] for i in pbin['discarded']], color='gray', s=2, zorder=5)
        else:
            self.ax.scatter(self.phase, self.denergy, color='magenta', s=4)
        plt.show()
        

    def update(self, num, redrawBackground=False):
        istart = num*self.nbr_p
        lost_particles = [i for i in self.active_particles if num in self.coll_hits and i in self.coll_hits[num]]
        self.active_particles = [i for i in self.active_particles if not i in lost_particles]

        # self.ref_e_text.set_text("E = {0:.8E}".format(self.ref_energy[num]))
        self.ref_e_text.set_text("Frame {}".format(num))
        self.pos_plot.set_offsets([(self.phase[i + istart], self.denergy[i + istart]) for i in self.active_particles])

        self.lost_particles['denergy'].extend((self.denergy[i + istart] for i in lost_particles))
        self.lost_particles['phase'].extend((self.phase[i + istart] for i in lost_particles))
        self.lost_plot.set_offsets([i for i in zip(self.lost_particles['phase'], self.lost_particles['denergy'])])

        if redrawBackground:
            try:
                filename = "phasespace/{}lines.dat".format(num)
                self.plot_trajectory(filename, redraw=True)
            except:
                print("could not redraw frame")

    def animate(self, save_to="", animateBackground=False):
        self.create_plot()
        self.plot_background_lines()
        self.plot_collimators()
        self.format_axes()

        self.lost_plot = self.ax.scatter([], [], color='r', marker='x', zorder=20, s=4)
        self.pos_plot = self.ax.scatter(self.phase[0:self.nbr_p], self.denergy[0:self.nbr_p], zorder=10, s=4, color='b')
        self.active_particles = range(self.nbr_p)
        self.lost_particles = {'denergy': [], 'phase': []}
        self.ref_e_text = self.ax.text(0.05, 0.95, 
                "Frame {}".format(0), fontsize=15,
                ha='left', va='top', transform=self.ax.transAxes)

        ani = animation.FuncAnimation(self.fig, self.update, self.nbr_turns, fargs=(animateBackground,), interval=50, blit=False)

        if save_to:
            print("saving simulation to '{}'".format(save_to))
            self.fig.suptitle(save_to)
            ani.save(save_to, fps=20, dpi=500)
            print("finished saving to '{}'".format(save_to))
        else:
            plt.show()

    def plot_energy_oscillations(self, particleID):
        raise Exception("deprecated -- removed support for 'self.ref_energy'. Read in 'ramp.txt' and convert it yourself if you want this.")

        turns = [i for i in range(self.nbr_turns) if i % 10 == 0]
        energy = [(e + self.denergy[i]) for i, e in enumerate(self.ref_energy) if i in turns]
        ref_e = [e for i, e in enumerate(self.ref_energy) if i in turns]

        fig, ax = plt.subplots()
        ax.plot(turns, energy, color='b')
        ax.plot(turns, ref_e, color='black')
        plt.show()


def plot_lossmap(lossmaps, labels=[], save_to=''):
    """
        lossmaps = [lossmap1, lossmap2, ...] 
    """
    if len(lossmaps) == 0:
        raise Exception("no losses found")

    ps = PhaseSpace(STARTDIST_FILE)
    max_turn = 0
    for lm in lossmaps:
        max_turn = max(max_turn, max(lm.keys()))
    turns = np.array(range(max_turn + 100))

    # Plotting
    fig, loss_ax = plt.subplots()

    # LOSSES
    # tot_loss = np.empty(len(turns))
    color_list = plt.cm.Set3(np.linspace(0, 1, len(lossmaps)))
    for i, lm in enumerate(lossmaps):
        # lm = lossmaps[action]
        losses = np.array([len(lm[turn]) if turn in lm else 0 for turn in turns])
        # tot_loss += losses
        avg_loss = moving_average(losses, int(1.3*11245))
        loss_ax.plot(turns, avg_loss, color=color_list[i], label=(labels[i] if len(labels) > i else ""), zorder=3)
    loss_ax.set_ylabel("Losses (∆particles/1.3s)")
    loss_ax.set_xlabel("t (s)")
    loss_ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.1f}".format(x/11245.0)))
    loss_ax.spines['right'].set_position(('axes', 1.15))
    loss_ax.set_yscale('log')
    loss_ax.legend(loc='upper right')

    # INTENSITY
    # n = ps.nbr_p
    # intensity = np.empty(len(tot_loss))
    # for i, v in enumerate(tot_loss):
        # n -= v
        # intensity[i] = float(n)/ps.nbr_p
    # intensity_ax = loss_ax.twinx()
    # intensity_ax.plot(turns, intensity, color='b', label='intensity', linestyle='--', zorder=1)
    # intensity_ax.set_ylabel("Intensity")

    # # RAMP
    ramp = np.array(read_ramp(RAMP_FILE, len(turns))['e'])
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

def export_fitted_ramp():
    ramp = read_ramp(RAMP_FILE, 500000)
    with open("ramp_fitted.txt", 'w') as f:
        for i, v in enumerate(ramp['e_fitted']):
            f.write("{} {:.16f}\n".format(i, v/1e6))

def export_particles(phasespace, plist, save_file):
    print("writing {} particles to '{}'".format(len(plist), save_file))
    with open(save_file, 'w') as f:
        f.write("{},1\n".format(len(plist)))
        for pid in plist:
            f.write("{},{},{}\n".format(phasespace.denergy[pid], phasespace.phase[pid], phasespace.h[pid]))


def plot_energy_oscillations():
    nbr_turns = 500*11245
    ramp = read_ramp(RAMP_FILE, nbr_turns)
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

def plot_hamiltonian(ps): # input phase space containing ps.h
    fig = plt.figure()
    ax = fig.add_subplot(111)

    lossmap = get_lossmap(COLL_FILE)
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

def plot_lost():
    # lossmap = get_lossmap(COLL_FILE, ["id"])
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

def plot_distribution(file_path):
    """ Plot the first frame of a given file. Will color particles depending if they eventually got lost or not """

    ps = PhaseSpace(file_path)
    # lossmap = get_lossmap(COLL_FILE, "id")
    # pbin = ps.categorize_particles(lossmap)
    ps.plot_particles()




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



if __name__ == "__main__":
    ACTION = argv[1]
    if ACTION == "animate":
        print("animate trajectory")
        ps = PhaseSpace(PARTICLE_FILE)
        ps.animate(save_to=SAVE_FILE)
    if ACTION == "animate-full":
        print("animating trajectory and background")
        output = "calc/animate-full.mp4"
        print("will save movie to '{}'".format(output))
        ps = PhaseSpace(PARTICLE_FILE)
        ps.animate(output, animateBackground=True)
    elif ACTION == "lossmap" or ACTION == "lossmap-analysis":
        print("plot lossmap")
        lossmap = get_lossmap(COLL_FILE)
        plot_lossmap([lossmap], SAVE_FILE)
    elif ACTION == "separated-lossmap":
        print("plot one series for each action value of the lossmap")
        ps = PhaseSpace(STARTDIST_FILE)
        lossmap = get_lossmap(COLL_FILE)
        div_lossmaps = {}
        for turn in lossmap:
            for pid in lossmap[turn]:
                h = int(ps.h[pid])
                if not h in div_lossmaps: 
                    div_lossmaps[h] = {}
                lossm = div_lossmaps[h]
                if not turn in lossm:
                    lossm[turn] = []
                lossm[turn].append(pid)
        # lossmaps = [div_lossmaps[action] for action in div_lossmaps]
        labels = sorted(div_lossmaps.keys())
        lossmaps = [div_lossmaps[l] for l in labels]
        plot_lossmap(lossmaps, labels)
    elif ACTION == "energy":
        print("plot energy oscillations")
        plot_energy_oscillations()
    elif ACTION == "startdist":
        print("plot start distribution")
        losses_from_sec = float(argv[2]) if len(argv) > 2 else 0.0
        ps = PhaseSpace(STARTDIST_FILE)
        lossmap = get_lossmap(COLL_FILE)
        pbin = ps.categorize_particles(lossmap, losses_from_sec)
        ps.plot_particles(pbin)
    elif ACTION == "dist":
        print("plot dist")
        if not len(argv) > 2: raise Exception("need to give a file path to plot")
        plot_distribution(argv[2])
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
    elif ACTION == "trajectory":
        print("plot particle trajectory")
        ps = PhaseSpace(PARTICLE_FILE)
        ps.create_plot()
        ps.format_axes()
        ps.plot_trajectory(randomizeColors=True)
        plt.show()
    elif ACTION == "ham":
        print("plot hamiltonian")
        ps = PhaseSpace(PARTICLE_FILE)
        plot_hamiltonian(ps)
    elif ACTION == "phasespace-mov":
        print("animating phasespace evolution")
        phasespace_evolution()
    elif ACTION == "lost":
        print("lost plot")
        plot_lost()
    else:
        print("unrecognised action")
