import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cmx
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle

import numpy as np
from settings import *

from math import pi

class PhaseSpace:
    """ Only works for data output by the 2d synchrotron
        
        Intended for two things:
            1. Conveniently store the data from the 2d synchrotron
            2. Plot things related to the phase space
        Thus, there is some mish-mash of functions below
    """

    @classmethod
    def merge_two(clss, ps1, ps2):
        """ ONLY WORKS FOR SINGLE DISTRIBUTIONS RIGHT NOW  (that is, two particles.dat files can't be merged)
            
            Can only merge two PhaseSpaces with the same number of turns.
        """
        if ps1 == None:
            return ps2
        elif ps2 == None:
            return ps1
        elif not ps1.nbr_turns == ps2.nbr_turns:
            raise Exception("both PhaseSpaces needs to have the same amount of turns")

        ps_merge = PhaseSpace(None)
        ps_merge.nbr_p = ps1.nbr_p + ps2.nbr_p
        ps_merge.nbr_turns = ps1.nbr_turns



        attributes = ['denergy', 'phase', 'h']
        for attr in attributes:
            a1 = getattr(ps1, attr)
            a2 = getattr(ps2, attr)
            setattr(ps_merge, attr, np.concatenate((a1, a2)))
            # a1 = getattr(ps1, attr).reshape(ps1.nbr_p, ps1.nbr_turns)
            # a2 = getattr(ps2, attr).reshape(ps2.nbr_p, ps2.nbr_turns)
            # setattr(ps_merge, attr, np.hstack((a1, a2)).reshape(ps_merge.nbr_p*ps_merge.nbr_turns))

        return ps_merge

    def __init__(self, pfile, mute=False):
        self.nbr_p = 0
        self.nbr_turns = 0
        self.denergy = np.array([])
        self.phase = np.array([])
        self.background_lines = []


        if pfile:
            if not mute: print("reading in particles from '{}'".format(pfile))
            with open(pfile, 'r') as file:
                self.nbr_p , self.nbr_turns = map(int, file.readline().strip().split(','))
                self.denergy = np.empty(self.nbr_p*self.nbr_turns)
                self.phase = np.empty(self.nbr_p*self.nbr_turns)
                self.h = np.empty(self.nbr_p*self.nbr_turns)
                for i, line in enumerate(file.readlines()):
                    denergy, phase, h = map(float, line.rstrip('\n').split(','))
                    self.denergy[i] = denergy
                    self.phase[i] = phase
                    self.h[i] = h

    def create_plot(self):
        self.fig, self.ax = plt.subplots()

    def plot_trajectory(self, filePath=None, randomizeColors=False, redraw=False):
        print("plot trajectory '{}'".format(filePath))
        phasespace = self
        if filePath is not None:
            phasespace = PhaseSpace(filePath)
        else:
            print("using local data")

        nbr_series = phasespace.nbr_p
        series = [{"denergy" : [], "phase" : [], "h" : 0.0} for l in range(nbr_series)]
        for i, de in enumerate(phasespace.denergy[:-1]):
            series[i % nbr_series]["phase"].append(phasespace.phase[i])
            series[i % nbr_series]["denergy"].append(de)
            series[i % nbr_series]['h'] = phasespace.h[i]

        if redraw and not len(self.background_lines) == len(series):
            raise Exception("new frame does not match series in old")

        cMap = plt.get_cmap('plasma_r')
        cNorm = colors.Normalize(vmin=0, vmax=4e6)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cMap)
        print("plotting series")
        color_list = plt.cm.Set3(np.linspace(0, 1, len(series)))
        for i, trail in enumerate(series):
            colorVal = scalarMap.to_rgba(trail['h'])
            if randomizeColors:
                colorVal = color_list[i]
            if redraw:
                self.background_lines[i].set_xdata(trail["phase"])
                self.background_lines[i].set_ydata(trail["denergy"])
            else:
                line, = self.ax.plot(trail["phase"], trail["denergy"], '-', color=colorVal, zorder=1)
                self.background_lines.append(line)

    def plot_background_lines(self):
        self.plot_trajectory(settings.LINE_PATH)

    def plot_collimators(self):
        print("plot collimators")
        try:
            with open(settings.COLL_PATH, 'r') as f:
                f.readline()
                second_line = f.readline()
                top_coll, bot_coll = map(float, second_line.rstrip().split(','))
                self.collimator = {'top' : top_coll, 'bot' : bot_coll}
        except:
            print("could not read '{}', will not plot".format(settings.COLL_PATH))
        else:
            # self.coll_hits = get_lossmap(COLL_FILE)
            self.ax.axhspan(settings.PLOT_FRAME['y'][0], self.collimator['bot'], facecolor='red', alpha=0.1)
            self.ax.axhspan(self.collimator['top'], settings.PLOT_FRAME['y'][1], facecolor='red', alpha=0.1)

    def format_axes(self):
        self.ax.set_xlim(settings.PLOT_FRAME['x'])
        self.ax.set_ylim(settings.PLOT_FRAME['y'])
        self.ax.set_xlabel("Phase (radians)")
        self.ax.set_ylabel("∆E (GeV)")
        self.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e9)))
        self.ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.1f}π".format(x/pi)))

    def categorize_particles(self, lossmap, sec=0.0):
        """ Seperates particle id's into three lists, 'lost' and 'alive' after time 'sec'.
            Particles are 'discarded' if they got lost before 'sec'. """
        pbin = {'lost' : [], 'discarded' : []}
        for turn in lossmap:
            if turn >= 11245.0*sec:
                pbin['lost'] += lossmap[turn]
            else:
                pbin['discarded'] += lossmap[turn]
        pbin['lost'] = np.array(sorted(pbin['lost']))
        pbin['discarded'] = np.array(sorted(pbin['discarded']))

        all_lost = np.concatenate((pbin['lost'], pbin['discarded']))
        all_p = np.arange(0, self.nbr_p, 1)

        pbin['alive'] = all_p[np.where(np.invert(np.in1d(all_p, all_lost)))]
        # The above is a _much_ quicker version of the (more readable) version below:
        # pbin['alive'] = [i for i in range(self.nbr_p) if i not in pbin['lost'] and i not in pbin['discarded']]
        return pbin

    def plot_particles(self, pbin=None):
        self.create_plot()
        self.plot_background_lines()
        self.plot_collimators()
        self.format_axes()
        
        if pbin:
            print("\talive:", len(pbin['alive']))
            print("\tlost:", len(pbin['lost']))
            print("\tdiscarded:", len(pbin['discarded']))
            print("\ttotal:", sum(map(len, pbin.values())))
            live_z = 10
            lost_z = 9 if len(pbin['lost']) > len(pbin['alive']) else 11
            if len(pbin['alive']) > 0: self.ax.scatter(self.phase[pbin['alive']], self.denergy[pbin['alive']], color='b', s=2, zorder=live_z, label='Alive')
            if len(pbin['lost']) > 0: self.ax.scatter(self.phase[pbin['lost']], self.denergy[pbin['lost']], color='r', s=2, zorder=lost_z, label='Lost')
            if len(pbin['discarded']) > 0: self.ax.scatter(self.phase[pbin['discarded']], self.denergy[pbin['discarded']], color='gray', s=2, zorder=5, label='Discarded')
            self.ax.legend(loc='lower left')
        else:
            self.ax.scatter(self.phase, self.denergy, s=4)
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

    def animate(self, lossmap, save_to="", animateBackground=False):
        self.coll_hits = lossmap

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
