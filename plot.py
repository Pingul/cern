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
PARTICLE_FILE = "data/particles.dat"
LINE_FILE = "data/lines.dat"
COLL_FILE = "data/coll.dat"
# RAMP_FILE = "resources/ramp.txt"
RAMP_FILE = "resources/cLHC_momentum_programme6.5TeV.dat"
STARTDIST_FILE = "data/startdist.dat"

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

# def get_lossmap(collfile, with_id=False, with_coll_values=False):
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

class Trajectory:
	def __init__(self):
		self.nbr_p = 0
		self.nbr_turns = 0
		self.denergy = []
		self.phase = []
		self.ref_energy = []

		contents = []
		print("reading in particles from '{}'".format(PARTICLE_FILE))
		with open(PARTICLE_FILE, 'r') as file:
			for i, line in enumerate(file.readlines()):
				if i == 0:
					self.nbr_p, self.nbr_turns = map(int, line.rstrip('\n').split(','))
					continue
				denergy, phase = map(float, line.rstrip('\n').split(','))
				self.denergy.append(denergy)
				self.phase.append(phase)

		print("reading in energy from '{}'".format(RAMP_FILE))
		with open(RAMP_FILE, 'r') as file:
			for i, line in enumerate(file.readlines()):
				if i >= self.nbr_turns:
					break
				ref_e = float(line.split()[-1])
				self.ref_energy.append(ref_e)

		self.fig, self.ax = plt.subplots()


	def plot_background_lines(self):
		nbr_lines = 0
		lines = []
		print("reading background lines from '{}'".format(LINE_FILE))
		with open(LINE_FILE, 'r') as f:
			for i, line in enumerate(f.readlines()):
				if i == 0:
					nbr_lines = int(line.rstrip().split(',')[0])
					lines = [{"denergy" : [], "phase" : []} for l in range(nbr_lines)]
					continue
				denergy, phase = map(float, line.rstrip("\n").split(","))
				lines[i % nbr_lines]["phase"].append(phase)
				lines[i % nbr_lines]["denergy"].append(denergy)

		cMap = plt.get_cmap('plasma_r')
		cNorm = colors.Normalize(vmin=0, vmax=6e8)
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cMap)
		print("plotting background lines")
		for trail in lines:
			mmax = abs(max(trail['denergy']))
			mmin = abs(min(trail['denergy']))
			v = mmax if mmax > mmin else mmin
			colorVal = scalarMap.to_rgba(v)
			self.ax.plot(trail["phase"], trail["denergy"], '-', color=colorVal, zorder=1)

	def fetch_collimators(self):
		with open(COLL_FILE, 'r') as f:
			f.readline()
			second_line = f.readline()
			top_coll, bot_coll = map(float, second_line.rstrip().split(','))
			self.collimator = {'top' : top_coll, 'bot' : bot_coll}

		coll_hits = get_lossmap(COLL_FILE, with_attr=['id'])
		self.coll_hits = {turn : [hit['id'] for hit in coll_hits[turn]] for turn in coll_hits}

		self.ax.axhspan(PLOT_FRAME['y'][0], self.collimator['bot'], facecolor='red', alpha=0.1)
		self.ax.axhspan(self.collimator['top'], PLOT_FRAME['y'][1], facecolor='red', alpha=0.1)

	def update(self, num):
		istart = num*self.nbr_p
		lost_particles = [i for i in self.active_particles if num in self.coll_hits and i in self.coll_hits[num]]
		self.active_particles = [i for i in self.active_particles if not i in lost_particles]

		self.ref_e_text.set_text("E = {0:.8E}".format(self.ref_energy[num]))
		self.pos_plot.set_offsets([(self.phase[i + istart], self.denergy[i + istart]) for i in self.active_particles])

		self.lost_particles['denergy'].extend((self.denergy[i + istart] for i in lost_particles))
		self.lost_particles['phase'].extend((self.phase[i + istart] for i in lost_particles))
		self.lost_plot.set_offsets([i for i in zip(self.lost_particles['phase'], self.lost_particles['denergy'])])

	def animate(self, save_to=""):
		self.ax.set_xlim(PLOT_FRAME['x'])
		self.ax.set_ylim(PLOT_FRAME['y'])
		self.ax.set_xlabel("Phase (radians)")
		self.ax.set_ylabel("∆E (GeV)")
		self.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:g}".format(x/1e9)))
		self.ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.1f}π".format(x/math.pi)))
		self.ref_e_text = self.ax.text(-4, 1.7e9, "E = {0:.4E}".format(self.ref_energy[0]), ha = 'left', va = 'center', fontsize = 15)

		self.lost_plot = self.ax.scatter([], [], color='r', marker='x', zorder=20)
		self.pos_plot = self.ax.scatter(self.phase[0:self.nbr_p], self.denergy[0:self.nbr_p], zorder=10)
		self.active_particles = range(self.nbr_p)
		self.lost_particles = {'denergy': [], 'phase': []}

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

def animate_trajectory():
	traj = Trajectory()
	traj.plot_background_lines()
	traj.fetch_collimators()
	traj.animate(save_to=SAVE_FILE)

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

def find_spikes(averaged_loss_data):
	cycles = []
	istart = -1
	iend = -1
	for i, d in enumerate(averaged_loss_data):
		if istart == -1 and d > 0:
			istart = i
		elif istart > 0 and d == 0:
			iend = i
			cycles.append([istart, iend])
			istart = iend = -1

	return cycles



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

	average_losses = moving_average(losses, 1125)

	# spikes = find_spikes(average_losses)
	print(len(spikes), "spikes:")
	for i, spike in enumerate(spikes, 1):
		print("\t{:g}: {:.2f} -> {:.2f} s ".format(i, secs[spike[0]], secs[spike[1]]))


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
	loss_ax.set_ylabel("Losses (∆particles/0.1s)")
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
	de = np.gradient(e)

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
	nbr_turns = 500000
	ramp = read_ramp(RAMP_FILE, nbr_turns)
	turns = range(nbr_turns)

	fig = plt.figure()
	e_ax = fig.add_subplot(211)
	de_ax = fig.add_subplot(212, sharex=e_ax)

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

	fig.suptitle("LHC ramp")

	plt.show()


if __name__ == "__main__":
	ACTION = argv[1]
	if ACTION == "animate":
		print("animate trajectory")
		animate_trajectory()
	elif ACTION == "lossmap" or ACTION == "lossmap-analysis":
		print("plot lossmap")
		plot_lossmap(SAVE_FILE)
		# plot_lossmap_phase()
	elif ACTION == "energy":
		print("plot energy oscillations")
		plot_energy_oscillations()
	else:
		print("unrecognised action")
