import os
import pytimber
import matplotlib.pyplot as plt
import math
import ast
import pickle
import numpy as np
import bisect

from scipy import stats
from scipy import interpolate

def store_file_for_fill(fill_nbr):
	return 'fills/fill_%s.dat' % str(fill_nbr)
	# return 'test_fill/fill_%s.dat' % str(fill_nbr)

def subset_indices(sequence, minv, maxv):
	low = bisect.bisect_left(sequence, minv)
	high = bisect.bisect_left(sequence, maxv, lo=low)
	return [low, high]

def moving_average(sequence, N):
	""" Moving average given the sequence. Returns an list equal in length
	to the one given """

	average = np.convolve(sequence, np.ones((N,))/N, mode='same')
	return average

def imax(data):
	""" Returns both the max and the index for the max value """
	i = max(range(len(data)), key = data.__getitem__)
	return (data[i], i)

class Fill:

	class STATUS:
		OK = 'OK'
		NORAMP = 'NORAMP'
		ERROR = 'ERROR'
	#### class STATUS

	class Variable:
		""" Purely a wrapper so we can get named variables instead of [0] and [1] """
		def __init__(self, data):
			# data = [[x1, x2, ..], [y1, y2, ...]]
			#
			# Using x and y because it's easy to write and quickly shows intent,
			# even though our x almost exclusively is time
			self.x = data[0]
			self.y = data[1]

		def __len__(self):
			return 2

		def __getitem__(self, key):
			if key == 0:
				return self.x
			elif key == 1:
				return self.y
			else:
				raise IndexError("not valid key '{}'".format(key))

		def __setitem__(self, key, value):
			if key == 0:
				self.time = x
			elif key == 1:
				self.val = y
			else:
				raise IndexError("not valid key '{}'".format(key))
	#### class Variable

	variables = {
		'intensity_b1' : 'LHC.BCTFR.A6R4.B1:BEAM_INTENSITY',
		# 'intensity_b2' : 'LHC.BCTFR.A6R4.B2:BEAM_INTENSITY',
		'energy' : 'LHC.BOFSU:OFSU_ENERGY',
		'beta_coll_b_b1' : 'BLMTI.06L7.B1E10_TCP.B6L7.B1:LOSS_RS09',
		'beta_coll_c_b1' : 'BLMTI.06L7.B1E10_TCP.C6L7.B1:LOSS_RS09',
		'beta_coll_d_b1' : 'BLMTI.06L7.B1E10_TCP.D6L7.B1:LOSS_RS09',
		'synch_coll_b1' : 'BLMTI.06L3.B1I10_TCP.6L3.B1:LOSS_RS09',
		'abort_gap_int_b1' : 'LHC.BSRA.US45.B1:ABORT_GAP_INTENSITY',
	}

	aligned_variables = {
		## ALIGNED DATA
		## Needs to have a A_ prefix
		'A_beta_coll_b_b1' : 'BLMTI.06L7.B1E10_TCP.B6L7.B1:LOSS_RS09',
		'A_beta_coll_c_b1' : 'BLMTI.06L7.B1E10_TCP.C6L7.B1:LOSS_RS09',
		'A_beta_coll_d_b1' : 'BLMTI.06L7.B1E10_TCP.D6L7.B1:LOSS_RS09',
		'A_synch_coll_b1' : 'BLMTI.06L3.B1I10_TCP.6L3.B1:LOSS_RS09',

	}
	all_variables = {**variables, **aligned_variables}

	db = pytimber.LoggingDB()

	def __init__(self, nbr, fetch=True):
		self.nbr = nbr
		self.data = {}
		self.meta = {}
		self.status = Fill.STATUS.OK
		for key in self.variables:
			self.data[key] = Fill.Variable([[], []])

		if fetch:
			try:
				self.fetch()
				# self.convert_to_using_Variable()
				self.normalize_intensity()
			except Exception as e:
				print("could not initialize fill {}:".format(self.nbr))
				print("\t", e)

	def fetch(self, forced=False, cache=True):
		store_file = store_file_for_fill(self.nbr)
		to_fetch = self.variables
		fetch_aligned = forced
		if not forced and os.path.isfile(store_file):
			print('loading {}'.format(self.nbr))
			with open(store_file, 'rb') as f:
				self.unpack(pickle.loads(f.read()))

			non_cached_variables = []
			for v in self.all_variables:
				if not v in self.data.keys():
					if v.startswith("A_"): # We have some unfetched aligned variable
						fetch_aligned = True
					else:
						non_cached_variables.append(v)
			if len(non_cached_variables) == 0 and not fetch_aligned:
				return
			# Fetching missing varaibles
			to_fetch = non_cached_variables

		## Forced fetch
		print('fetching {}'.format(self.nbr))

		self.meta = self.db.getLHCFillData(self.nbr)
		start_t = self.meta['startTime']
		end_t = self.meta['endTime']
		preramp = next((item for item in self.meta['beamModes'] if item['mode'] == 'PRERAMP'), None)
		ramp = next((item for item in self.meta['beamModes'] if item['mode'] == 'RAMP'), None)
		if not preramp or not ramp:
			self.status = Fill.STATUS.NORAMP
		else:
			start_t = preramp['startTime']
			end_t = ramp['endTime']

		if to_fetch:
			print("\tfetch", to_fetch.keys())
			data = self.db.get(list(to_fetch.values()), start_t, end_t)
			for d_var in data:
				var_name = next(name for name in self.variables.keys() if self.variables[name] == d_var)
				self.data[var_name] = Fill.Variable(data[d_var])

		if fetch_aligned:
			# If we need to fetch a new aligned variable, we need to fetch all old
			# ones as well for all to be aligned
			print("\tfetch aligned", self.aligned_variables.keys())
			data = self.db.getAligned(list(self.aligned_variables.values()), start_t, end_t)
			x_data = data.pop('timestamps')
			for d_var in data:
				var_name = next(name for name in self.aligned_variables.keys() if self.aligned_variables[name] == d_var)
				self.data[var_name] = Fill.Variable(np.array([x_data, data[d_var]]))

		if cache:
			print('caching {}'.format(self.nbr))
			with open(store_file, 'wb') as f:
				pickle.dump(self.pack(), f)

	def pack(self):
		return {
			'nbr' : self.nbr,
			'data' : self.data,
			'meta' : self.meta,
			'status' : self.status
		}

	def unpack(self, dump):
		self.nbr = dump['nbr']
		self.data = dump['data']
		self.meta = dump['meta']
		self.status = dump['status']

	def normalize_intensity(self):
		for v in ['intensity_b1']:
			m = max(self.data[v].y)
			if m > 0.0:
				self.data[v].y = self.data[v].y/m

	def has_off_momentum_loss(self, beam='b1'):
		variable = 'synch_coll_%s' % beam
		off_momentum_losses = np.array(self.data[variable].y)
		# d_off_momentum_losses = np.gradient(off_momentum_losses)

		max_loss = max(off_momentum_losses)
		mean_loss = np.mean(off_momentum_losses)
		if max_loss < 0 or mean_loss < 0:
			return False
		return max_loss/mean_loss > 10

	def beta_coll_merge(self):
		# Should really be dealing with aligned data here
		beta_var = ['A_beta_coll_b_b1', 'A_beta_coll_c_b1', 'A_beta_coll_d_b1', ]
		if not False in [i in self.data.keys() for i in beta_var]:
			beta_coll = None
			for v in beta_var:
				if beta_coll is None:
					beta_coll = self.data[v]
				else:
					beta_coll.y = np.add(beta_coll.y, self.data[v].y)
			self.data['A_beta_coll_b1'] = beta_coll



	## PLOTTING
	def plot(self):
		print('plotting {}'.format(self.nbr))
		fig = plt.figure()
		intensity_axis  = fig.add_subplot(211)
		energy_axis = intensity_axis.twinx()
		blm_axis = intensity_axis.twinx()


		intensity_axis.plot(*self.data['intensity_b1'], color='b', zorder=10, linestyle='-', linewidth='1')
		intensity_axis.set_ylim([0.95, 1.005])
		intensity_axis.set_ylabel("Beam Intensity")

		energy_axis.plot(*self.data['energy'], color='black', zorder=5)
		energy_axis.set_ylabel("Energy")

		oml = self.data['synch_coll_b1']
		spike, tail = find_spike(oml.y)
		blm_axis.axvspan(oml.x[spike], oml.x[tail], facecolor='b', alpha=0.2)

		fig.subplots_adjust(right=0.8)
		blm_axis.spines['right'].set_position(('axes', 1.15))
		blm_axis.set_frame_on(True)
		blm_axis.plot(*oml, color='r', linestyle='--', zorder=1)
		blm_axis.set_yscale('log')
		blm_axis.set_ylabel("Losses")


		plt.title("Fill {}".format(self.nbr))

		agap_ax = fig.add_subplot(212, sharex=intensity_axis)
		agap_ax.plot(self.data['abort_gap_int_b1'].x, moving_average(self.data['abort_gap_int_b1'].y, 10), color='g')
		agap_ax.set_ylabel("Abort gap intensity")

		plt.show()

	def blm_plot(self):
		print('loss plot {}'.format(self.nbr))
		self.beta_coll_merge()

		fig = plt.figure()
		intensity_axis = fig.add_subplot(111)
		energy_axis = intensity_axis.twinx()
		blm_axis = intensity_axis.twinx()

		intensity_axis.plot(*self.data['intensity_b1'], color='b', zorder=10, linestyle='-', linewidth='1')
		intensity_axis.set_ylim([0.95, 1.005])
		intensity_axis.set_ylabel("Beam Intensity")

		energy_axis.plot(*self.data['energy'], color='black', zorder=5)
		energy_axis.set_ylabel("Energy")

		oml = self.data['synch_coll_b1']
		spike, tail = find_spike(oml.y)

		fig.subplots_adjust(right=0.8)
		blm_axis.spines['right'].set_position(('axes', 1.15))
		blm_axis.set_frame_on(True)
		blm_axis.plot(*oml, color='r', linestyle='--', zorder=2, label='momentum loss')
		blm_axis.plot(*self.data['A_beta_coll_b1'], color='g', linestyle='--', zorder=1, label='transversal loss')
		blm_axis.axvspan(oml.x[spike], oml.x[tail], facecolor='b', alpha=0.2)
		blm_axis.set_yscale('log')
		blm_axis.set_ylabel("Losses")
		blm_axis.legend(loc='lower right')

		plt.title("Fill {}".format(self.nbr))
		plt.show()

	def oml_plot(self):
		print('loss plot {}'.format(self.nbr))
		self.beta_coll_merge()

		fig = plt.figure()
		intensity_axis = fig.add_subplot(111)
		energy_axis = intensity_axis.twinx()
		blm_axis = intensity_axis.twinx()

		intensity_axis.plot(*self.data['intensity_b1'], color='b', zorder=10, linestyle='-', linewidth='1')
		intensity_axis.set_ylim([0.95, 1.005])
		intensity_axis.set_ylabel("Beam Intensity")

		energy_axis.plot(*self.data['energy'], color='black', zorder=5)
		energy_axis.set_ylabel("Energy")

		fig.subplots_adjust(right=0.8)
		blm_axis.spines['right'].set_position(('axes', 1.15))
		blm_axis.set_frame_on(True)

		oml = self.data['synch_coll_b1']
		avg = moving_average(oml.y, 10)
		spike, tail = find_spike(oml.y)
		blm_axis.plot(*oml, color='g', linestyle='--', zorder=2, label='momentum loss')
		blm_axis.plot(oml.x, avg, color='r', linestyle='-', zorder=2, label='average momentum loss')
		blm_axis.axvspan(oml.x[spike], oml.x[tail], facecolor='b', alpha=0.2)

		blm_axis.set_yscale('log')
		blm_axis.set_ylabel("Losses")
		blm_axis.legend(loc='lower right')

		plt.title("Fill {}".format(self.nbr))
		plt.show()

	def convert_to_using_Variable(self):
		store_file = store_file_for_fill(self.nbr)
		new_data = {}
		for v in self.data.keys():
			new_data[v] = Fill.Variable(self.data[v])
		self.data = new_data
		with open(store_file, 'wb') as f:
			pickle.dump(self.pack(), f)


	### DEPRECATED
	def crop_ramp(self):
		energy = np.array(self.data['energy'][1])
		d_energy = np.gradient(energy)

		index_limit = {'max' : 0, 'min' : len(d_energy)}
		# Sometimes a fill contains multiple ramps. We want the last one,
		# and use the extra trial index limit to help us find it
		trial_index_limit = index_limit.copy()

		time_limit = {'max' : 0, 'min' : 0}
		threshold = 4 # seems to work to find the ramp
		for index, de in enumerate(d_energy):
			if de > threshold and index < index_limit['min']:
				index_limit['min'] = index
			elif de < threshold and index > index_limit['min']:
				index_limit['max'] = index
				trial_index_limit = index_limit.copy()
				index_limit = {'max' : 0, 'min' : len(d_energy)}

		if not trial_index_limit['max'] > trial_index_limit['min']:
			raise Exception("could not find ramp")
		index_limit = trial_index_limit

		padding = 500
		time_limit['max'] = self.data['energy'][0][index_limit['max']] + padding
		time_limit['min'] = self.data['energy'][0][index_limit['min']] - padding

		for key in self.data:
			low, high = subset_indices(self.data[key][0], time_limit['min'], time_limit['max'])
			self.data[key] = [l[low:high] for l in self.data[key]]

##### class Fill
################





def evaluate_off_momentum_losses_in_fills(fills, save_file):
	open(save_file, 'w').close() # erasing file

	for i in fills:
		print("evaluating %s..." % str(i))
		fill = Fill(i, fetch=False)
		fill.fetch(forced=False)
		with open(save_file, 'a') as f:
			f.write(str(i))
			status_string = '\t'
			if fill.status == Fill.STATUS.OK:
				status_string += 'OML' if fill.has_off_momentum_loss() else 'OK'
			else:
				status_string += str(fill.status)
			status_string += '\n'
			f.write(status_string)
		print("--\n")

def fills_from_file(file, status_string='*'):
	fills = []
	for line in open(file, 'r'):
		contents = line.rstrip('\n').split()
		fill_nbr = int(contents[0])
		status = contents[1]

		if status_string == "*" or status.upper() == status_string.upper():
			fills.append(fill_nbr)
	print("Found {} fills".format(len(fills)))
	return fills

def plot_from(file, status_string='*'):
	fills = fills_from_file(file, status_string)
	plot_at_the_time = 10
	n = 0
	for fill_nbr in fills:
		n += 1
		print("evaluating %s" % fill_nbr)
		fill = Fill(fill_nbr, fetch=False)
		fill.fetch(forced=False, cache=False)
		fill.normalize_intensity()
		# if not status_string.upper() == 'ERROR':
		# 	fill.crop_ramp()
		fill.plot()
		print("--\n")
		if n % plot_at_the_time == 0:
			inp = input("draw {} more plots? (press 'q' to quit) ".format(plot_at_the_time))
			if inp == 'q':
				break

def plot_energy_ramp(fill_nbr):
	fill = Fill(fill_nbr)

	energy = fill.data['energy']

	fig = plt.figure()
	ax1 = fig.add_subplot(211)
	plt.title("Energy ramp and derivative around start of ramp for fill 5427")

	ax1.plot(*energy)
	ax1.set_ylabel("Energy GeV")

	ax2 = fig.add_subplot(212, sharex=ax1)
	d_energy = np.gradient(energy.y)/np.gradient(energy.x)
	ax2.plot(energy.x, d_energy)
	ax2.set_ylabel("∆ Energy GeV")

	plt.show()

def find_start_of_ramp(energy): 
	denergy = np.gradient(energy.y)
	dt = np.gradient(energy.x)
	denergyDt = denergy/dt # the spacing could be non-uniform, so we need this to correctly represent the derivative

	istart = 0
	e_threshold = 450.0
	de_threshold = 0.1
	for i, e in enumerate(energy.y):
		de = denergyDt[i]
		if e > e_threshold and de > de_threshold:
			istart = i
			break
	else:
		raise Exception("could not find the start of the ramp")
	return istart


def find_spike(data): # data is normally BLM data from momentum collimators
	ddata = np.abs(np.gradient(data))
	mov_average = moving_average(ddata, 10)
	threshold = 1e-7

	vpeak, ipeak = imax(data)
	start = end = ipeak
	found_start = found_end = False
	while not found_start and start > 0:
		start -= 1
		if ddata[start] < threshold and data[start] < vpeak/5e1:
			found_start = True

	while not found_end and end < len(ddata) - 1:
		end += 1
		if ddata[end] < threshold and data[end] < vpeak/1e2:
			found_end = True

	# if not found_start or not found_end:
	# 	raise Exception("could not find ({})spike")

	return [start, end]


def find_crossover_point(betaBLM, momentumBLM):
	""" Look for the point after OML spike when transversal losses starts 
		to dominate the momentum losses 

		Note: this should be used with aligned data """
	vpeak, ipeak = imax(momentumBLM.y)

	i = ipeak
	while betaBLM[i].y < momentumBLM[i].y:
		i += 1
		if i >= len(betaBLM.y):
			raise Exception("could not find crossover point")
	return i



def draw_histogram(title, data, binsize, xlabel='', ylabel='', color='b'):
	maxbin = max(data) + binsize
	minbin = min(data)
	bins = np.arange(minbin, maxbin, binsize)
	fig, ax = plt.subplots()
	ax.hist(data, bins=bins, color=color)
	ax.set_title(title)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	plt.show()

def intensity_and_OML_pruning(file_in, file_out):
	fills = fills_from_file(file_in, "OML")

	open(file_out, 'w').close() # erasing file

	low_intensity = 0
	wrongly_categorised = 0
	for nbr in fills:
		fill = Fill(nbr, False)
		fill.fetch()

		if max(fill.data['intensity_b1'].y) < 1.8e14:
			low_intensity += 1
			continue

		fill.beta_coll_merge()
		oml_data = fill.data['synch_coll_b1']
		smin, smax = find_spike(oml_data.y) 
		dur = oml_data.x[smax] - oml_data.x[smin]

		if dur < 100 or dur > 300:
			wrongly_categorised += 1
			continue

		with open(file_out, 'a') as f:
			f.write("{}\tOML\n".format(str(fill.nbr)))

	removed = int(round(float(low_intensity + wrongly_categorised)/len(fills) * 100))
	print('discarded:')
	print('\tintensity {}'.format(low_intensity))
	print('\tmis-categorised {}'.format(wrongly_categorised))
	print('\ttotal {}%'.format(removed))

def export_energy_ramp_to_sixtrack(file_out='ramp.txt', fill=5433):
	f = Fill(fill)
	freq = 11245 # laps / s
	with open(file_out, 'w') as file:
		x = f.data['energy'].x[0:300] # at least 300 s, should be enough
		y = f.data['energy'].y[0:300]
		rescaled_x = [freq*(i - x[0]) for i in x]
		delta = int(round(rescaled_x[-1] - rescaled_x[0]))

		new_x = np.linspace(0, delta, delta)
		new_y = interpolate.interp1d(rescaled_x, y, kind='cubic')(new_x)
		for i, e in enumerate(new_y):
			scaled_e = e*1e3 # Sixtrack uses the unit MeV
			turnnbr = int(round(new_x[i])) + 1 # 1 indexed
			file.write('{} {}\n'.format(turnnbr, scaled_e))


#################
## Statistics

def intensity_histogram(file):
	fills = fills_from_file(file, "OML")
	intensities = []
	for nbr in fills:
		fill = Fill(nbr, fetch=False)
		fill.fetch()
		# fill.crop_ramp()

		intensities.append(max(fill.data['intensity_b1'].y))

	draw_histogram('Intensity for {}'.format(file), intensities, 1e13, "Intensity", "Count")

def spike_duration_histogram(file):
	fills = fills_from_file(file, "OML")
	outliers = []
	durations = []
	for nbr in fills:
		fill = Fill(nbr)

		losses = fill.data['synch_coll_b1'].y
		spike_start, spike_end = find_spike(np.array(losses))
		d = fill.data['synch_coll_b1'].x[spike_end] - fill.data['synch_coll_b1'].x[spike_start]
		if d < 70 or d > 300:
			outliers.append(nbr)
		durations.append(d)

	draw_histogram('Spike duration for {}'.format(file), durations, 10, 'Seconds', 'Count')
	return outliers

def spike_energy_histogram(file):
	fills = fills_from_file(file, "OML")
	spike_energy = []
	spike_time = []

	for nbr in fills:
		fill = Fill(nbr)

		losses = fill.data['synch_coll_b1']
		ispike = np.argmax(losses.y)
		tspike = losses.x[ispike]

		energy = fill.data['energy']
		ienergy = bisect.bisect_left(energy.x, tspike)

		iramp = find_start_of_ramp(energy)
		delta_e = energy.y[ienergy] - energy.y[iramp]
		delta_t = energy.x[ienergy] - energy.x[iramp]
		spike_energy.append(delta_e)
		spike_time.append(delta_t)

	draw_histogram('Spike energy for {}'.format(file), spike_energy, 0.13, 'Delta E (GeV) from start of ramp', 'Count', 'y')
	draw_histogram('Max spike event for {}'.format(file), spike_time, 1.0, 'Delta t (s) from start of ramp till spike', 'Count')

def beta_vs_synch_blm(file):
	fills = fills_from_file(file, "OML")
	ok = 0
	notok = 0
	sdata = {
		'max' : [],
		'mean' : [],
	}
	bdata = {
		'max' : [],
		'mean' : []
	}
	for nbr in fills:
		fill = Fill(nbr)
		try:
			fill.beta_coll_merge()
			ok +=1
		except Exception as e:
			notok += 1
			continue

		smin, smax = find_spike(fill.data['synch_coll_b1'].y) 
		tmax = fill.data['synch_coll_b1'].x[smax]
		tmin = fill.data['synch_coll_b1'].x[smin]
		bmin, bmax = subset_indices(fill.data['A_beta_coll_b1'].x, tmin, tmax)

		bsubset = fill.data['A_beta_coll_b1'].y[bmin:bmax]
		ssubset = fill.data['synch_coll_b1'].y[smin:smax]

		sdata['max'].append(max(ssubset))
		sdata['mean'].append(np.mean(ssubset))
		bdata['max'].append(max(bsubset))
		bdata['mean'].append(np.mean(bsubset))

	fig, ax = plt.subplots()
	ax.set_xlabel("Synchrotron (IR3) TCP")
	ax.set_ylabel("Betatron (IR7) TCPs")

	ax.scatter(sdata['max'], bdata['max'], color='r')
	slope, intercept, r_value, p_value, std_err = stats.linregress(sdata['max'], bdata['max'])
	print(slope, intercept, r_value, p_value, std_err)
	xval = [0, 1]
	max_yval = [slope*x + intercept for x in xval]
	ax.plot(xval, max_yval, color='r', label='max')

	ax.scatter(sdata['mean'], bdata['mean'], color='b')
	slope, intercept, r_value, p_value, std_err = stats.linregress(sdata['mean'], bdata['mean'])
	print(slope, intercept, r_value, p_value, std_err)
	mean_yval = [slope*x + intercept for x in xval]
	ax.plot(xval, mean_yval, color='b', label='mean')

	ax.plot([0, 1], [0, 1], color = 'black', label='delimiter')

	for v in ['max', 'mean']:
		count = 0
		for i, sd in enumerate(sdata[v]):
			if bdata[v][i] > sd: 
				count += 1
		print(v, "over: ", count, "({}%)".format(int(float(count)/len(sdata[v])*100)))

	files_used = int(float(ok)/(ok + notok) * 100)
	plt.title('Losses due to synchrotron vs betatron oscillations\n for {} (could use {}% of the fills)'.format(file, files_used))
	ax.legend(loc='upper right')
	ax.set_ylim([0, 0.5])
	ax.set_xlim([0, 0.5])
	plt.show()

def intensity_vs_OML(file):
	fills = fills_from_file(file, "OML")
	intensity = []
	mean_loss = []
	max_loss = []
	discarded = 0
	for nbr in fills:
		fill = Fill(nbr, False)
		fill.fetch()
		smin, smax = find_spike(fill.data['synch_coll_b1'].y) 
		ssubset = fill.data['synch_coll_b1'].y[smin:smax]

		maxint = max(fill.data['intensity_b1'][1])
		if maxint < 1.8e14:
			discarded += 1
			continue

		mean_loss.append(np.mean(ssubset))
		max_loss.append(max(ssubset))
		intensity.append(maxint)

	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122, sharey=ax1) 

	# fig1, ax1 = plt.subplots()
	ax1.set_xlabel("Mean momentum (IR3) TCP")
	ax1.set_ylabel("Intensity")
	ax1.scatter(mean_loss, intensity, color='b', label='mean')
	ax1.set_xlim([0, 1.1*max(mean_loss)])
	ax1.set_ylim([1.5e14, 1.1*max(intensity)])
	ax1.legend(loc="lower right")

	# fig2, ax2 = plt.subplots()
	ax2.set_xlabel("Max momentum (IR3) TCP")
	ax2.set_ylabel("Intensity")
	ax2.scatter(max_loss, intensity, color='r', label='max')
	ax2.set_xlim([0, 1.1*max(max_loss)])
	ax2.legend(loc="lower right")

	percent_used = int(round(float(len(intensity))/(len(intensity) + discarded) * 100))
	fig.suptitle("Intensity vs OML for {} (only intenities > 1.8e14, {}% of total)\n".format(file, percent_used))

	plt.show()

def abort_gap_vs_OML(file):
	fills = fills_from_file(file, "OML")
	abort_gap = []
	mean_loss = []
	max_loss = []
	for nbr in fills:
		fill = Fill(nbr, False)
		fill.fetch()
		smin, smax = find_spike(fill.data['synch_coll_b1'].y) 
		tmax = fill.data['synch_coll_b1'].x[smax]
		tmin = fill.data['synch_coll_b1'].x[smin]
		agmin, agmax = subset_indices(fill.data['abort_gap_int_b1'].x, tmin, tmax)
		ag_average = moving_average(fill.data['abort_gap_int_b1'].y, 5)

		ssubset = fill.data['synch_coll_b1'].y[smin:smax]

		mean_loss.append(np.mean(ssubset))
		max_loss.append(max(ssubset))
		abort_gap.append(ag_average[agmin] - ag_average[agmax])


	fig = plt.figure()
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122, sharey=ax1) 

	# fig1, ax1 = plt.subplots()
	ax1.set_xlabel("Mean momentum (IR3) TCP")
	ax1.set_ylabel("Abort gap intensity")
	ax1.scatter(mean_loss, abort_gap, color='b', label='mean')
	ax1.set_xlim([0, 1.1*max(mean_loss)])
	ax1.set_ylim([0, 1.1*max(abort_gap)])
	ax1.legend(loc="lower right")

	# fig2, ax2 = plt.subplots()
	ax2.set_xlabel("Max momentum (IR3) TCP")
	ax2.set_ylabel("Abort gap intensity")
	ax2.scatter(max_loss, abort_gap, color='r', label='max')
	ax2.set_xlim([0, 1.1*max(max_loss)])
	ax2.legend(loc="lower right")

	fig.suptitle("Abort gap intensity vs OML for {}\n".format(file))

	plt.show()
	
def abort_gap_vs_BLM(file):
	fills = fills_from_file(file, "OML")
	abort_gap = []
	smean_loss = []
	bmean_loss = []
	for nbr in fills:
		fill = Fill(nbr, False)
		fill.fetch()
		fill.beta_coll_merge()

		smin, smax = find_spike(fill.data['synch_coll_b1'].y) 
		tmax = fill.data['synch_coll_b1'].x[smax]
		tmin = fill.data['synch_coll_b1'].x[smin]
		agmin, agmax = subset_indices(fill.data['abort_gap_int_b1'].x, tmin, tmax)
		ag_average = moving_average(fill.data['abort_gap_int_b1'].y, 5)

		ssubset = fill.data['synch_coll_b1'].y[smin:smax]
		bsubset = fill.data['A_beta_coll_b1'].y[smin:smax]

		smean_loss.append(np.mean(ssubset))
		bmean_loss.append(np.mean(bsubset))
		abort_gap.append(ag_average[agmin] - ag_average[agmax])

	tot_mean_loss = np.add(np.array(smean_loss), np.array(bmean_loss))
	vmax = 1.1*max(tot_mean_loss)

	fig = plt.figure()
	ax1 = fig.add_subplot(131)
	ax2 = fig.add_subplot(132, sharey=ax1) 
	ax3 = fig.add_subplot(133, sharey=ax1) 

	# fig1, ax1 = plt.subplots()
	ax1.set_xlabel("Momentum loss")
	ax1.set_ylabel("Abort gap intensity")
	ax1.scatter(smean_loss, abort_gap, color='r', label='longitudinal')
	ax1.set_xlim([0, vmax])
	ax1.xaxis.set_ticks(np.arange(0, vmax, 0.005))
	ax1.set_ylim([0, 1.1*max(abort_gap)])
	# ax1.legend(loc="lower right")

	ax2.set_xlabel("Transversal loss")
	ax2.scatter(bmean_loss, abort_gap, color='g', label='transversal')
	ax2.set_xlim([0, vmax])
	ax2.xaxis.set_ticks(np.arange(0, vmax, 0.005))

	ax3.set_xlabel("Total loss")
	ax3.scatter(tot_mean_loss, abort_gap, color='b', label='total')
	ax3.set_xlim([0, vmax])
	ax3.xaxis.set_ticks(np.arange(0, vmax, 0.005))

	fig.suptitle("Abort gap intensity vs BLM for {}\n".format(file))

	plt.show()
		


## I think these are irrelevant. Up for deletion.
##
# def loss_int_histogram(file):
# 	fills = fills_from_file(file, "OML")
# 	integrated_losses = []
# 	for nbr in fills:
# 		fill = Fill(nbr)

# 		losses = np.array(fill.data['synch_coll_b1'][1])
# 		spike, tail = find_spike(losses)
# 		int_loss = 0.0
# 		for i in range(spike, tail):
# 			int_loss += losses[i]
# 		integrated_losses.append(int_loss)

# 	draw_histogram('Integrated losses', integrated_losses, 0.01)
# 	return integrated_losses


# def max_loss_histogram(file):
# 	fills = fills_from_file(file, "OML")
# 	max_loss = []
# 	for nbr in fills:
# 		fill = Fill(nbr, fetch=False)
# 		fill.fetch()
# 		# fill.crop_ramp() # can throw
# 		max_loss.append(max(fill.data['synch_coll_b1'][1]))

# 	draw_histogram('Max losses', max_loss, 0.005)
# 	return max_loss

# def find_spike(data):
# 	raise Exception("should not be used")
# 	# ddata = np.gradient(data)
# 	spike_index = max(range(len(data)), key=data.__getitem__) # - 3 # just so we get more of the start of the spike
# 	# spike_val = data[spike_index]

# 	# Used 21 previously --> feels unreasonable
# 	mov_average = moving_average(data, 5)
# 	# dmov_average = np.gradient(mov_average)

# 	# Use the moving average to calculate the tail
# 	average_max = max(mov_average)
# 	tail_threshold = average_max/10.0
# 	tail_index = -1
# 	for i, d in enumerate(mov_average[spike_index:]):
# 		if d < tail_threshold:
# 			tail_index = spike_index + i
# 			break
# 	else:
# 		print("Spike: ", spike_index, tail_index)
# 		raise Exception("did not find tail")

# 	# Note that
# 	#    spike_index < tail_index
# 	return [spike_index, tail_index]

