import os
import pytimber
import matplotlib.pyplot as plt
import math
import ast
import pickle
import numpy as np
import bisect

def store_file_for_fill(fill_nbr):
	return 'fills/fill_%s.dat' % str(fill_nbr)

def subset_indices(sequence, minv, maxv):
	low = bisect.bisect_left(sequence, minv)
	high = bisect.bisect_left(sequence, maxv, lo=low)
	return [low, high]

def moving_average(sequence, N):
	""" Moving average given the sequence. Returns an list equal in length
	to the one given """

	average = np.convolve(sequence, np.ones((N,))/N, mode='same')
	return average

class Fill:
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
	db = pytimber.LoggingDB()

	def __init__(self, nbr, fetch=True):
		self.nbr = nbr
		self.data = {}
		self.meta = {}
		for key in self.variables:
			self.data[key] = []

		if fetch:
			try:
				self.fetch()
				self.normalize_intensity()
				self.crop_ramp()
			except Exception as e:
				print("could not initialize fill {}:".format(self.nbr))
				print("\t", e)

	def fetch(self, forced=False, cache=True):
		store_file = store_file_for_fill(self.nbr)

		to_fetch = self.variables
		if not forced and os.path.isfile(store_file):
			print('loading {}'.format(self.nbr))
			with open(store_file, 'rb') as f:
				self.data = pickle.loads(f.read())

			non_cached_variables = []
			for v in self.variables:
				if not v in self.data.keys():
					non_cached_variables.append(v)

			if len(non_cached_variables) == 0:
				return

			# Fetching missing varaibles
			to_fetch = non_cached_variables

		## Forced fetch
		print('fetching {}'.format(self.nbr))

		self.meta = self.db.getLHCFillData(self.nbr)
		start_t = self.meta['startTime']
		end_t = self.meta['endTime']

		for vkey in to_fetch:
			print('\tfetching ' + vkey)
			d = self.db.get(self.variables[vkey], start_t, end_t)
			for dkey in d:
				self.data[vkey] = list(d[dkey])

		# beta_variables = ['beta_coll_b_b1', 'beta_coll_c_b1', 'beta_coll_d_b1']
		# if not False in [i in self.data.keys for i in beta_variables]:
		# 	# Merging all betatron losses
		# 	self.data['beta_coll_b1'] = 
		# 	for v in beta_variables:


		if cache:
			print('caching {}'.format(self.nbr))
			with open(store_file, 'wb') as f:
				pickle.dump(self.data, f)


	def normalize_intensity(self):
		# for v in ['intensity_b1', 'intensity_b2']:
		for v in ['intensity_b1']:
			m = max(self.data[v][1])
			if m > 0.0:
				self.data[v][1] = [i/m for i in self.data[v][1]]

	def crop_ramp(self):
		# This would maybe have been a smarter method, but now it's just more work to re-fetch all data
		#
		# preramp = next((item for item in self.meta['beamModes'] if item['mode'] == 'PRERAMP'), None)
		# ramp = next((item for item in self.meta['beamModes'] if item['mode'] == 'RAMP'), None)
		# if not preramp or not ramp:
		# 	raise Exception("could not find ramp")
		# ramp_start = preramp['startTime']
		# ramp_end = ramp['endTime']

		# for key in self.data:
		# 	low, high = subset_indices(self.data[key][0], ramp_start, ramp_end)
		# 	self.data[key] = [l[low:high] for l in self.data[key]]

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

		fig.subplots_adjust(right=0.8)
		blm_axis.spines['right'].set_position(('axes', 1.15))
		blm_axis.set_frame_on(True)
		blm_axis.plot(*self.data['synch_coll_b1'], color='r', linestyle='--', zorder=1)
		blm_axis.set_yscale('log')
		blm_axis.set_ylabel("Losses")

		plt.title("Fill {}".format(self.nbr))

		agap_ax = fig.add_subplot(212, sharex=intensity_axis)
		agap_ax.plot(self.data['abort_gap_int_b1'][0], moving_average(self.data['abort_gap_int_b1'][1], 10), color='g')
		agap_ax.set_ylabel("Abort gap intensity")

		plt.show()

	def blm_plot(self):
		print('loss plot {}'.format(self.nbr))
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
		blm_axis.plot(*self.data['synch_coll_b1'], color='r', linestyle='--', zorder=2)
		blm_axis.plot(*self.data['beta_coll_b1'], color='g', linestyle='--', zorder=1)
		blm_axis.set_yscale('log')
		blm_axis.set_ylabel("Losses")

		plt.title("Fill {}".format(self.nbr))
		plt.show()


	def has_off_momentum_loss(self, beam='b1'):
		variable = 'synch_coll_%s' % beam
		off_momentum_losses = np.array(self.data[variable][1])
		d_off_momentum_losses = np.gradient(off_momentum_losses)

		max_loss = max(d_off_momentum_losses)
		mean_loss = np.mean(d_off_momentum_losses)
		if max_loss < 0 or mean_loss < 0:
			return False
		return max_loss/mean_loss > 1000

	def beta_coll_merge(self):
		beta_var = ['beta_coll_b_b1', 'beta_coll_c_b1', 'beta_coll_d_b1', ]
		if not False in [i in self.data.keys() for i in beta_var]:
			beta_coll = self.data[beta_var[0]]
			for i in range(2, len(beta_var)):
				beta_coll[1] = np.add(beta_coll[1], self.data[beta_var[i]][1])
			self.data['beta_coll_b1'] = beta_coll


# END OF FILL

def evaluate_off_momentum_losses_in_fills(fills, save_file):
	open(save_file, 'w').close() # erasing file

	for i in fills:
		print("evaluating %s..." % str(i))
		fill = Fill(i, fetch=False)
		fill.fetch(forced=False)
		with open(save_file, 'a') as f:
			f.write(str(i))
			status_string = '\t'
			try:
				fill.crop_ramp()
				status_string += 'OML' if fill.has_off_momentum_loss() else 'OK'
			except Exception as e:
				status_string += 'NORAMP'
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
		if not status_string.upper() == 'ERROR':
			fill.crop_ramp()
		fill.plot()
		print("--\n")
		if n % plot_at_the_time == 0:
			inp = input("draw {} more plots? (press 'q' to quit) ".format(plot_at_the_time))
			if inp == 'q':
				break


def find_spike(data):
	ddata = np.gradient(data)
	spike_index = max(range(len(ddata)), key=ddata.__getitem__) # - 3 # just so we get more of the start of the spike
	spike_val = data[spike_index]

	# moving average on data
	padding = 10
	start = spike_index
	end = len(data) - padding
	moving_average = []
	for i in range(start, end):
		mean = 0
		for j in range(i - padding, i + padding):
			mean += data[j]
		mean /= 2.0*float(padding)
		moving_average.append(mean)

	# Use the moving average to calculate the tail
	average_max = max(moving_average)
	tail_threshold = average_max/10.0
	tail_index = -1
	for i, d in enumerate(moving_average):
		if d < tail_threshold:
			tail_index = spike_index + i
			break
	else:
		print("Spike: ", spike_index, tail_index)
		raise Exception("did not find tail")

	return [spike_index, tail_index]

def test():
	fill = Fill(5433)
	# losses = np.array(fill.data['synch_coll_b1'])
	# spike, tail = find_spike(losses[1])
	# print(spike, tail)
	# fig, ax = plt.subplots()
	fig, blm_axis = plt.subplots()
	ax2 = blm_axis.twinx()

	ax2.plot(*fill.data['abort_gap_int_b1'], zorder=1)

	# fig.subplots_adjust(right=0.8)
	# blm_axis.spines['right'].set_position(('axes', 1.15))
	blm_axis.set_frame_on(True)
	blm_axis.plot(*fill.data['synch_coll_b1'], color='r', linestyle='--', zorder=10)
	# blm_axis.plot(*self.data['beta_coll_b1'], color='b', linestyle='-.')
	blm_axis.set_yscale('log')
	blm_axis.set_ylabel("Losses")


	plt.show()


def draw_histogram(title, data, binsize):
	maxbin = max(data) + binsize
	bins = np.arange(0, maxbin, binsize)
	fig, ax = plt.subplots()
	ax.hist(data, bins=bins)
	ax.set_title(title)
	plt.show()




#################
## Statistics

def acc_loss_vs_intensity(file):
	fills = fills_from_file(file, "OML")
	intensities = []
	integrated_losses = []
	for nbr in fills:
		fill = Fill(nbr, fetch=False)
		fill.fetch()
		fill.crop_ramp()

		intensities.append(max(fill.data['intensity_b1'][1]))

		losses = np.array(fill.data['synch_coll_b1'][1])
		spike, tail = find_spike(losses)
		int_loss = 0.0
		for i in range(spike, tail):
			int_loss += losses[i]
		integrated_losses.append(int_loss)

	fig, ax = plt.subplots()
	ax.scatter(intensities, integrated_losses)
	ax.set_ylabel("Integrated loss")
	ax.set_xlabel("Intensity")
	plt.title("Loss vs intensity for {}".format(file))
	plt.show()

def intensity_histogram(file):
	fills = fills_from_file(file, "OML")
	intensities = []
	for nbr in fills:
		fill = Fill(nbr, fetch=False)
		fill.fetch()
		fill.crop_ramp()

		intensities.append(max(fill.data['intensity_b1'][1]))

	draw_histogram('Intensity', intensities, 1e13)

def loss_duration_histogram(file):
	fills = fills_from_file(file, "OML")
	durations = []
	for nbr in fills:
		fill = Fill(nbr)

		losses = fill.data['synch_coll_b1'][1]
		spike, tail = find_spike(np.array(losses))
		d = fill.data['synch_coll_b1'][0][tail] - fill.data['synch_coll_b1'][0][spike]
		durations.append(d)

	draw_histogram('Duration of spike', durations, 5)

def loss_int_histogram(file):
	fills = fills_from_file(file, "OML")
	integrated_losses = []
	for nbr in fills:
		fill = Fill(nbr)

		losses = np.array(fill.data['synch_coll_b1'][1])
		spike, tail = find_spike(losses)
		int_loss = 0.0
		for i in range(spike, tail):
			int_loss += losses[i]
		integrated_losses.append(int_loss)

	draw_histogram('Integrated losses', integrated_losses, 0.01)
	return integrated_losses


def max_loss_histogram(file):
	fills = fills_from_file(file, "OML")
	max_loss = []
	for nbr in fills:
		fill = Fill(nbr, fetch=False)
		fill.fetch()
		fill.crop_ramp() # can throw
		max_loss.append(max(fill.data['synch_coll_b1'][1]))

	draw_histogram('Max losses', max_loss, 0.005)
	return max_loss

