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
		print(max_loss, mean_loss, max_loss/mean_loss)
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

		fig.subplots_adjust(right=0.8)
		blm_axis.spines['right'].set_position(('axes', 1.15))
		blm_axis.set_frame_on(True)
		blm_axis.plot(*self.data['synch_coll_b1'], color='r', linestyle='--', zorder=1)
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

def intensity_and_OML_pruning(file_in, file_out):
	fills = fills_from_file(file_in, "OML")

	open(file_out, 'w').close() # erasing file

	for nbr in fills:
		fill = Fill(nbr, False)
		fill.fetch()

		if max(fill.data['intensity_b1'].y) < 1.8e14:
			continue

		fill.beta_coll_merge()
		smin, smax = find_spike(fill.data['synch_coll_b1'].y) 
		tmax = fill.data['synch_coll_b1'].x[smax]
		tmin = fill.data['synch_coll_b1'].x[smin]
		bmin, bmax = subset_indices(fill.data['A_beta_coll_b1'].x, tmin, tmax)

		bsubset = fill.data['A_beta_coll_b1'].y[bmin:bmax]
		ssubset = fill.data['synch_coll_b1'].y[smin:smax]

		if np.mean(bsubset) > np.mean(ssubset):
			continue

		with open(file_out, 'a') as f:
			f.write("{}\tOML\n".format(str(fill.nbr)))





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

		intensities.append(max(fill.data['intensity_b1'].y))

		losses = np.array(fill.data['synch_coll_b1'].y)
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

		intensities.append(max(fill.data['intensity_b1'].y))

	draw_histogram('Intensity for {}'.format(file), intensities, 1e13, "Intensity", "Count")

def loss_duration_histogram(file):
	fills = fills_from_file(file, "OML")
	durations = []
	for nbr in fills:
		fill = Fill(nbr)

		losses = fill.data['synch_coll_b1'].y
		spike, tail = find_spike(np.array(losses))
		d = fill.data['synch_coll_b1'].x[tail] - fill.data['synch_coll_b1'].x[spike]
		durations.append(d)

	draw_histogram('Spike duration for {}'.format(file), durations, 5, 'Seconds', 'Count')

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

		delta_e = energy.y[ienergy] - min(energy.y)
		spike_energy.append(delta_e)

		ramp = next((item for item in fill.meta['beamModes'] if item['mode'] == 'RAMP'), None)
		# Sometimes data is weirdly labeled, so we take the latest time of the below
		start_t = max(ramp['startTime'], energy.x[0])
		delta_t = energy.x[ienergy] - start_t
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

