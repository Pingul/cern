import os
import pytimber
import matplotlib.pyplot as plt
import math
import ast
import pickle
import numpy as np
import bisect

from scipy import stats

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

	class STATUS:
		OK = 'OK'
		NORAMP = 'NORAMP'
		ERROR = 'ERROR'
	#### class STATUS

	# class Variable:
	# 	def __init__(self, data):
	# 		# data = [[t1, t2, ..], [v1, v2, ...]]
	# 		self.time

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
		self.status = Fill.STATUS.OK
		for key in self.variables:
			self.data[key] = []

		if fetch:
			try:
				self.fetch()
				self.normalize_intensity()
				# self.crop_ramp()
			except Exception as e:
				print("could not initialize fill {}:".format(self.nbr))
				print("\t", e)

	def fetch(self, forced=False, cache=True):
		self.__fetch__('get', forced, cache)

	def fetchAligned(self):
		# Because we might lose data due to this, we fetch it newly every time.
		# Could possibly try to cache it seperately
		self.__fetch__('getAligned', forced=True, cache=False)

	def __fetch__(self, func, forced=False, cache=True):
		store_file = store_file_for_fill(self.nbr)
		to_fetch = self.variables
		if not forced and os.path.isfile(store_file):
			print('loading {}'.format(self.nbr))
			with open(store_file, 'rb') as f:
				self.unpack(pickle.loads(f.read()))
				# self.data = pickle.loads(f.read())

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
		preramp = next((item for item in self.meta['beamModes'] if item['mode'] == 'PRERAMP'), None)
		ramp = next((item for item in self.meta['beamModes'] if item['mode'] == 'RAMP'), None)
		if not preramp or not ramp:
			self.status = Fill.STATUS.NORAMP
		else:
			start_t = preramp['startTime']
			end_t = ramp['endTime']

		for vkey in to_fetch:
			print('\tfetching ' + vkey)
			func_call = 'self.db.{}(self.variables["{}"], {}, {})'.format(func, vkey, start_t, end_t)
			d = eval(func_call)
			# d = self.db.func(self.variables[vkey], start_t, end_t)
			for dkey in d:
				self.data[vkey] = list(d[dkey])

		if cache:
			print('caching {}'.format(self.nbr))
			with open(store_file, 'wb') as f:
				pickle.dump(self.pack(), f)
				# pickle.dump(self.data, f)

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
		# for v in ['intensity_b1', 'intensity_b2']:
		for v in ['intensity_b1']:
			m = max(self.data[v][1])
			if m > 0.0:
				self.data[v][1] = self.data[v][1]/m
				# self.data[v][1] = [i/m for i in self.data[v][1]]

	def has_off_momentum_loss(self, beam='b1'):
		variable = 'synch_coll_%s' % beam
		off_momentum_losses = np.array(self.data[variable][1])
		# d_off_momentum_losses = np.gradient(off_momentum_losses)

		max_loss = max(off_momentum_losses)
		mean_loss = np.mean(off_momentum_losses)
		print(max_loss, mean_loss, max_loss/mean_loss)
		if max_loss < 0 or mean_loss < 0:
			return False
		return max_loss/mean_loss > 10

	def beta_coll_merge(self):
		beta_var = ['beta_coll_b_b1', 'beta_coll_c_b1', 'beta_coll_d_b1', ]
		if not False in [i in self.data.keys() for i in beta_var]:
			# Looking for the largest possible subsequence
			# tmin = max([self.data[i][0][0] for i in beta_var])
			# tmax = min([self.data[i][0][-1] for i in beta_var])

			# print(tmin, tmax)

			beta_coll = None
			for v in beta_var:
				if beta_coll is None:
					beta_coll = self.data[v]
				else:
					beta_coll[1] = np.add(beta_coll[1], self.data[v][1])
			self.data['beta_coll_b1'] = beta_coll



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
		blm_axis.plot(*self.data['synch_coll_b1'], color='r', linestyle='--', zorder=2, label='momentum loss')
		blm_axis.plot(*self.data['beta_coll_b1'], color='g', linestyle='--', zorder=1, label='transversal loss')
		blm_axis.set_yscale('log')
		blm_axis.set_ylabel("Losses")
		blm_axis.legend(loc='lower right')

		plt.title("Fill {}".format(self.nbr))
		plt.show()



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

	# Note that
	#    spike_index < tail_index
	return [spike_index, tail_index]

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




#################
## Statistics

def acc_loss_vs_intensity(file):
	fills = fills_from_file(file, "OML")
	intensities = []
	integrated_losses = []
	for nbr in fills:
		fill = Fill(nbr, fetch=False)
		fill.fetch()
		# fill.crop_ramp()

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
		# fill.crop_ramp()

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

	draw_histogram('Spike duration for {}'.format(file), durations, 5, 'Seconds', 'Count')

def spike_energy_histogram(file):
	fills = fills_from_file(file, "OML")
	spike_energy = []
	spike_time = []

	for nbr in fills:
		fill = Fill(nbr)

		losses = fill.data['synch_coll_b1']
		ispike = np.argmax(losses[1])
		tspike = losses[0][ispike]

		energy = fill.data['energy']
		ienergy = bisect.bisect_left(energy[0], tspike)

		delta_e = energy[1][ienergy] - min(energy[1])
		spike_energy.append(delta_e)

		ramp = next((item for item in fill.meta['beamModes'] if item['mode'] == 'RAMP'), None)
		# Sometimes data is weirdly labeled, so we take the latest time of the below
		start_t = max(ramp['startTime'], energy[0][0])
		delta_t = energy[0][ienergy] - start_t
		spike_time.append(delta_t)


	draw_histogram('Spike energy for {}'.format(file), spike_energy, 0.1, 'Delta E (GeV) from start of ramp', 'Count', 'y')
	draw_histogram('Max spike event for {}'.format(file), spike_time, 2, 'Delta t (s) from start of ramp till spike', 'Count')

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

		smin, smax = find_spike(fill.data['synch_coll_b1'][1]) 
		tmax = fill.data['synch_coll_b1'][0][smax]
		tmin = fill.data['synch_coll_b1'][0][smin]
		# print(fill.data['beta_coll_b1'][0], tmin, tmax)
		bmin, bmax = subset_indices(fill.data['beta_coll_b1'][0], tmin, tmax)

		bsubset = fill.data['beta_coll_b1'][1][bmin:bmax]
		ssubset = fill.data['synch_coll_b1'][1][smin:smax]

		sdata['max'].append(max(ssubset))
		sdata['mean'].append(np.mean(ssubset))
		bdata['max'].append(max(bsubset))
		bdata['mean'].append(np.mean(bsubset))

	for v in ['max', 'mean']:
		sdata[v].sort()
		bdata[v].sort()

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

	files_used = int(float(ok)/(ok + notok) * 100)
	plt.title('Losses due to synchrotron vs betatron oscillations\n for {} (could use {}% of the fills)'.format(file, files_used))
	ax.legend(loc='upper right')
	ax.set_ylim([0, 0.5])
	ax.set_xlim([0, 0.5])
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

