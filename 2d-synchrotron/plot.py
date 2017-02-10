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
    'x' : [-2*math.pi, 4*math.pi],
    'y' : [-2e9, 2e9]
}

SAVE_FILE = ''
if (len(argv) > 2):
    SAVE_FILE = argv[2]

def moving_average(sequence, N):
    """ Moving average given the sequence. Returns an list equal in length
    to the one given """

    average = np.convolve(sequence, np.ones((N,))/N, mode='same')
    return average

def get_lossmap(collfile, with_attr=['id', 'phase', 'e']):
    # the default attributes are the only one that gets recognised
    # no values only gives nbr of losses each turn

    print("creating lossmap")
    hits = []
    with open(collfile, 'r') as f:
        print("reading collimation file '{}'".format(collfile))
        for i, line in enumerate(f):
            if i < 2: 
                continue
            line_c = line.rstrip().split(',')
            pID, turn = map(int, line_c[0:2])
            phase, e = map(float, line_c[2:4])
            hits.append([turn, pID, phase, e])

    hits.sort(key=lambda x: x[0])
    coll_hits = {} 
    for hit in hits:
        turn, pID, phase, e = hit

        val = {}
        if 'id' in with_attr:
            val['id'] = pID
        if 'phase' in with_attr:
            val['phase'] = phase
        if 'e' in with_attr:
            val['e'] = e

        try:
            coll_hits[turn].append(val)
        except Exception as e:
            coll_hits[turn] = [val]

    if not with_attr:
        return {turn : len(coll_hits[turn]) for turn in coll_hits}
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

    def plot_trajectory(self, filePath=None, randomizeColors=False):
        print("plot trajectory '{}'".format(filePath))
        nbr_series = 0
        series = []
        if filePath is not None:
            print("reading file '{}'".format(filePath))
            with open(filePath, 'r') as f:
                for i, line in enumerate(f.readlines()):
                    if i == 0:
                        nbr_series = int(line.rstrip().split(',')[0])
                        series = [{"denergy" : [], "phase" : []} for l in range(nbr_series)]
                        continue
                    denergy, phase, h = map(float, line.rstrip("\n").split(","))
                    series[i % nbr_series]["phase"].append(phase)
                    series[i % nbr_series]["denergy"].append(denergy)
        else:
            print("using local data")
            nbr_series = self.nbr_p
            series = [{"denergy" : [], "phase" : []} for l in range(nbr_series)]
            for i, de in enumerate(self.denergy):
                series[i % nbr_series]["phase"].append(self.phase[i])
                series[i % nbr_series]["denergy"].append(de)


        cMap = plt.get_cmap('plasma_r')
        cNorm = colors.Normalize(vmin=0, vmax=6e8)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cMap)
        print("plotting series")
        for trail in series:
            mmax = abs(max(trail['denergy']))
            mmin = abs(min(trail['denergy']))
            v = mmax if mmax > mmin else mmin
            colorVal = scalarMap.to_rgba(v)
            if randomizeColors:
                colorVal = np.random.rand(3,)
            self.ax.plot(trail["phase"], trail["denergy"], '-', color=colorVal, zorder=1)

    def plot_background_lines(self):
        self.plot_trajectory(LINE_FILE)

    def plot_collimators(self):
        with open(COLL_FILE, 'r') as f:
            f.readline()
            second_line = f.readline()
            top_coll, bot_coll = map(float, second_line.rstrip().split(','))
            self.collimator = {'top' : top_coll, 'bot' : bot_coll}

        coll_hits = get_lossmap(COLL_FILE, with_attr=['id'])
        self.coll_hits = {turn : [hit['id'] for hit in coll_hits[turn]] for turn in coll_hits}

        self.ax.axhspan(PLOT_FRAME['y'][0], self.collimator['bot'], facecolor='red', alpha=0.1)
        self.ax.axhspan(self.collimator['top'], PLOT_FRAME['y'][1], facecolor='red', alpha=0.1)

    def format_axes(self):
        self.ax.set_xlim(PLOT_FRAME['x'])
        self.ax.set_ylim(PLOT_FRAME['y'])
        self.ax.set_xlabel("Phase (radians)")
        self.ax.set_ylabel("∆E (GeV)")
        self.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e9)))
        self.ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.1f}π".format(x/math.pi)))

    def plot_turn(self, turn=1):
        if turn > self.nbr_turns or turn < 0:
            raise Exception("can't plot given turn")

        self.create_plot()
        self.plot_background_lines()
        self.plot_collimators()
        self.format_axes()

        lossmap = get_lossmap(COLL_FILE, "id")
        if not lossmap:
            start = (turn - 1)*self.nbr_p
            end = turn*self.nbr_p
            self.pos_plot = self.ax.scatter(self.phase[start:end], self.denergy[start:end], zorder=10, color='b', s=4)
        else:
            lost_particles = []
            for turn in lossmap:
                lost_particles += [loss['id'] for loss in lossmap[turn]]
            lost_particles.sort()
            live_particles = [i for i in range(self.nbr_p) if i not in lost_particles]

            livez = 10
            lostz = 9 if len(lost_particles) > len(live_particles) else 11
            self.ax.scatter([self.phase[i] for i in live_particles], [self.denergy[i] for i in live_particles], color='b', s=4, zorder=livez)
            self.ax.scatter([self.phase[i] for i in lost_particles], [self.denergy[i] for i in lost_particles], color='r', s=4, zorder=lostz)
        plt.show()

    def update(self, num):
        istart = num*self.nbr_p
        lost_particles = [i for i in self.active_particles if num in self.coll_hits and i in self.coll_hits[num]]
        self.active_particles = [i for i in self.active_particles if not i in lost_particles]

        # self.ref_e_text.set_text("E = {0:.8E}".format(self.ref_energy[num]))
        self.ref_e_text.set_text("Frame {}".format(num))
        self.pos_plot.set_offsets([(self.phase[i + istart], self.denergy[i + istart]) for i in self.active_particles])

        self.lost_particles['denergy'].extend((self.denergy[i + istart] for i in lost_particles))
        self.lost_particles['phase'].extend((self.phase[i + istart] for i in lost_particles))
        self.lost_plot.set_offsets([i for i in zip(self.lost_particles['phase'], self.lost_particles['denergy'])])

    def animate(self, save_to=""):
        self.create_plot()
        self.plot_background_lines()
        self.plot_collimators()
        self.format_axes()

        self.lost_plot = self.ax.scatter([], [], color='r', marker='x', zorder=20, s=4)
        self.pos_plot = self.ax.scatter(self.phase[0:self.nbr_p], self.denergy[0:self.nbr_p], zorder=10, s=4, color='b')
        self.active_particles = range(self.nbr_p)
        self.lost_particles = {'denergy': [], 'phase': []}
        # self.ref_e_text = self.ax.text(-4, 1.7e9, "E = {0:.4E}".format(self.ref_energy[0]), ha = 'left', va = 'center', fontsize = 15)
        self.ref_e_text = self.ax.text(-4, 1.7e9, "Frame {}".format(0), ha = 'left', va = 'center', fontsize = 15)

        ani = animation.FuncAnimation(self.fig, self.update, int(len(self.denergy)/self.nbr_p), interval=50, blit=False)

        if save_to:
            print("saving simulation to '{}'".format(save_to))
            self.fig.suptitle(save_to)
            ani.save(save_to, fps=20)
            print("finished saving")
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

def plot_lossmap_phase():
    lossmap = get_lossmap(COLL_FILE, with_attr=['id', 'phase'])

    turns = []
    phase = []
    pid = []

    for turn in lossmap:
        for loss in lossmap[turn]:
            turns.append(turn)
            phase.append(loss['phase'])
            pid.append(loss['id'])

    fig, ax = plt.subplots()
    ax.scatter(turns, phase)
    ax.set_xlabel("Time (kturns)")
    ax.set_ylabel("Phase")
    ax.set_xlim([0, max(lossmap.keys()) + 100])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1000.0)))
    fig.suptitle("Lossmap time vs phase")

    plt.show()


def plot_lossmap(save_to=''):
    lossmap = get_lossmap(COLL_FILE, with_attr=[])

    if len(lossmap) == 0:
        raise Exception("no losses found")

    max_turn_loss = max(lossmap.keys())
    nbr_p = 0
    coll = {'top' : 0, 'bot' : 0}
    with open(COLL_FILE, 'r') as f:
        first_line = f.readline()
        nbr_p = int(first_line.rstrip())
        second_line = f.readline()
        coll['top'], coll['bot'] = map(float, second_line.rstrip().split(','))


    turns = np.array(range(max_turn_loss + 100))
    secs = np.array([turn/11245.0 for turn in turns])
    losses = np.array([lossmap[turn] if turn in lossmap else 0 for turn in turns])
    ramp = np.array(read_ramp(RAMP_FILE, len(turns))['e'])

    n = nbr_p
    intensity = np.empty(len(losses))
    for i, v in enumerate(losses):
        n -= v
        intensity[i] = float(n)/nbr_p

    average_losses = moving_average(losses, int(1.3*11245))
    # average_losses = moving_average(losses, 1125)


    # Plotting
    fig, intensity_ax = plt.subplots()
    fig.subplots_adjust(right=0.8)

    intensity_ax.plot(turns, intensity, color='b', label='intensity')
    intensity_ax.set_xlabel("Time (kturns)")
    intensity_ax.set_ylabel("Intensity")
    intensity_ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    # intensity_ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.1f}".format(x/11245.0)))
    intensity_ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1000.0)))

    e_ax = intensity_ax.twinx()
    ramp_line = e_ax.plot(turns, ramp, color='black', zorder=0, label='LHC energy ramp')[0]
    e_ax.set_ylabel("E (GeV)")
    e_ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e9)))
    top_coll = Rectangle((0, 0), 1, 1, fc='w', fill=False, edgecolor='none', linewidth=0)
    bot_coll = Rectangle((0, 0), 1, 1, fc='w', fill=False, edgecolor='none', linewidth=0)
    e_ax.legend([ramp_line, top_coll, bot_coll], (ramp_line.get_label(), "Top coll:∆{:.2E}".format(coll['top']), "Bot coll:∆{:.2E}".format(coll['bot'])), loc='upper right')

    loss_ax = intensity_ax.twinx()
    loss_ax.plot(turns, average_losses, color='r', linestyle='--', label='∆loss')
    loss_ax.set_ylabel("Losses (∆particles/1.3s)")
    loss_ax.spines['right'].set_position(('axes', 1.15))
    if max(losses) > 0:
        loss_ax.set_yscale('log')

    fig.suptitle("Toy model lossmap")
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

    #e_ax.axvline(x=50*11245, linewidth=1)
    #de_ax.axvline(x=50*11245, linewidth=1)
    #ph_ax.axvline(x=50*11245, linewidth=1)

    e_ax.axvline(x=150*11245, linewidth=1, color='r')
    de_ax.axvline(x=150*11245, linewidth=1, color='r')
    ph_ax.axvline(x=150*11245, linewidth=1, color='r')

    fig.suptitle("LHC ramp")

    plt.show()

def plot_hamiltonian(ps): # input phase space containing ps.h
    fig = plt.figure()
    ax = fig.add_subplot(111)

    print(ps.h.shape)
    x = np.empty(ps.h.shape)
    for i in range(ps.nbr_turns):
        for j in range(ps.nbr_p):
            x[i*j + j] = i

    print("drawing {} points...".format(len(x)))
    ax.scatter(x, ps.h, color='brown')
    plt.show()



if __name__ == "__main__":
    ACTION = argv[1]
    if ACTION == "animate":
        print("animate trajectory")
        ps = PhaseSpace(PARTICLE_FILE)
        ps.animate(save_to=SAVE_FILE)
    elif ACTION == "lossmap" or ACTION == "lossmap-analysis":
        print("plot lossmap")
        plot_lossmap(SAVE_FILE)
        # plot_lossmap_phase()
    elif ACTION == "energy":
        print("plot energy oscillations")
        plot_energy_oscillations()
    elif ACTION == "startdist":
        print("plot start distribution")
        ps = PhaseSpace(STARTDIST_FILE)
        ps.plot_turn()
    elif ACTION == "enddist":
        print("plot end distribution")
        ps = PhaseSpace(ENDDIST_FILE)
        ps.plot_turn()
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
    else:
        print("unrecognised action")
