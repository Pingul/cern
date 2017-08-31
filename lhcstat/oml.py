import os, sys
import pytimber
import numpy as np
import pickle
import json

import bisect
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy import stats, interpolate, optimize

sys.path.append("../common")
from logger import ModuleLogger, LogLevel
from settings import settings

lg = ModuleLogger("oml")

FILLS = [5433, 5427, 5426, 5424, 5423, 5421, 5418, 5416, 5406, 5405, 5401, 5395, 5394, 5393, 5391, 5370, 5355, 5351, 5345, 5340, 5339, 5288, 5287, 5282, 5279, 5276, 5275, 5274, 5270, 5267, 5266, 5265, 5264, 5261, 5258, 5257, 5256, 5254, 5253, 5251, 5229, 5222, 5219, 5211, 5210, 5209, 5206, 5205, 5199, 5198, 5197, 5196, 5187, 5183, 5181, 5179, 5170, 5169, 5163, 5162, 5161, 5154, 5117, 5116, 5112, 5111, 5109, 5108, 5107, 5106, 5105, 5102, 5101, 5097, 5096, 5095, 5093, 5091, 5085, 5083, 5080, 5078, 5076, 5073, 5072, 5071, 5069, 5068, 5060, 5059, 5056, 5052, 5048, 5045, 5043, 5038, 5030, 5029, 5028, 5027, 5026, 5024, 5021, 5020, 5017, 5013, 4990, 4988, 4985, 4984, 4980, 4979, 4976, 4965, 4964, 4961, 4960]

def store_file_for_fill(fill_nbr):
    return os.path.join(os.path.dirname(__file__), 'fills/fill_{}.dat'.format(fill_nbr))

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
                self.x = value
            elif key == 1:
                self.y = value
            else:
                raise IndexError("not valid key '{}'".format(key))

        def index_for_time(self, timestamps):
            """ Helper function: find the index in the variable corresponding to the given timestamp """
            try:
                indices = []
                for t in timestamps:
                    indices.append(self.index_for_time(t))
                return indices
            except TypeError:
                index = bisect.bisect_left(self.x, timestamps)

                # Sometimes there's jump in the data sets -- take the value closest to the given timestamp
                ts = np.abs(self.x[[index-1, index]] - timestamps)
                return index-1 if ts[0] < ts[1] else index
    #### class Variable

    ## Variables in lists will be fetched aligned with each other
    BEAM1_VARIABLES = {
        'intensity_b1' : 'LHC.BCTFR.A6R4.B1:BEAM_INTENSITY',
        # 'beta_coll_b1' : ['BLMTI.06L7.B1E10_TCP.B6L7.B1:LOSS_RS09', 'BLMTI.06L7.B1E10_TCP.C6L7.B1:LOSS_RS09', 'BLMTI.06L7.B1E10_TCP.D6L7.B1:LOSS_RS09'],
        'beta_coll_b1' : 'BLMTI.06L7.B1E10_TCP.C6L7.B1:LOSS_RS09', 
        'synch_coll_b1' : 'BLMTI.06L3.B1I10_TCP.6L3.B1:LOSS_RS09',
        'abort_gap_int_b1' : 'LHC.BSRA.US45.B1:ABORT_GAP_INTENSITY',

        'motor_ir3_b1' : ['TCP.6L3.B1:MEAS_MOTOR_LU', 'TCP.6L3.B1:MEAS_MOTOR_RU'],
        'motor_ir7_b1' : ['TCP.C6L7.B1:MEAS_MOTOR_LU', 'TCP.C6L7.B1:MEAS_MOTOR_RU'],
    }

    BEAM2_VARIABLES = {
        'intensity_b2' : 'LHC.BCTFR.A6R4.B2:BEAM_INTENSITY',  
        # 'beta_coll_b2' : ['BLMTI.06R7.B2I10_TCP.B6R7.B2:LOSS_RS09', 'BLMTI.06R7.B2I10_TCP.C6R7.B2:LOSS_RS09', 'BLMTI.06R7.B2I10_TCP.D6R7.B2:LOSS_RS09'],
        'beta_coll_b2' : 'BLMTI.06R7.B2I10_TCP.C6R7.B2:LOSS_RS09', 
        'synch_coll_b2' : 'BLMTI.06R3.B2E10_TCP.6R3.B2:LOSS_RS09', 
        'abort_gap_int_b2' : 'LHC.BSRA.US45.B2:ABORT_GAP_INTENSITY',

        'motor_ir3_b2' : ['TCP.6R3.B2:MEAS_MOTOR_LU', 'TCP.6R3.B2:MEAS_MOTOR_RU'],
        'motor_ir7_b2' : ['TCP.C6R7.B2:MEAS_MOTOR_LU', 'TCP.C6R7.B2:MEAS_MOTOR_RU'],
    }

    COMB_VARIABLES = {
        'energy' : 'LHC.BOFSU:OFSU_ENERGY',
        'ramp_mode' : 'HX:BMODE_RAMP',
        'motor_start' : 'TCP.C6L7.B1:MEAS_PROFILE_TIME',
    }
    

    db = pytimber.LoggingDB()

    OML_period_file = "fills/spikes.dat"

    def __init__(self, nbr, fetch=True, beam=settings.BEAM):
        """ Note -- the default beam is decided upon loading this file """

        if not beam in (1, 2):
            raise Exception("Beam can only be 1 or 2")
        self.nbr = nbr
        self.beam = beam
        self.data = {}
        self.meta = {}
        self.status = Fill.STATUS.OK
        self.time_correction = 0

        if self.beam == 1:
            self.timber_var_map = {**self.BEAM1_VARIABLES, **self.COMB_VARIABLES}
        else:
            self.timber_var_map = {**self.BEAM2_VARIABLES, **self.COMB_VARIABLES}

        if fetch:
            self.fetch()
            self.normalize_intensity()
            self.beta_coll_merge()
            self.offset_time()


    def fetch(self, forced=False, cache=True):
        self.fetch_range("PRERAMP", "RAMP", forced, cache)

    def fetch_range(self, start, stop, forced=False, cache=True):
        l = ["INJPROT", "INJPHYS", "PRERAMP", "RAMP", "FLATTOP", "SQUEEZE", "ADJUST", "STABLE", "BEAMDUMP", "RAMPDOWN"]
        if not start in l or not stop in l:
            raise Exception("fetch range [{} - {}] is not allowed".format(start, stop))

        to_fetch = self.timber_var_map.keys()
        cache_file = store_file_for_fill(self.nbr)
        if not forced and os.path.isfile(cache_file):
            self.load_cache()
            cached_variables = list(self.data.keys())
            non_cached_variables = [var for var in self.timber_var_map.keys() if not var in cached_variables]
            to_fetch = non_cached_variables

        if not to_fetch:
            return

        aligned_var = [var for var in to_fetch if type(self.timber_var_map[var]) == list]
        non_aligned_var = [var for var in to_fetch if type(self.timber_var_map[var]) == str]
        if not (len(aligned_var) + len(non_aligned_var)) == len(to_fetch): 
            # Sanity check
            raise Exception("Sanity check failed")

        self.meta = self.db.getLHCFillData(self.nbr)
        start_t = self.meta['startTime']
        end_t = self.meta['endTime']
        start_mode = next((item for item in self.meta['beamModes'] if item['mode'] == start), None)
        end_mode = next((item for item in self.meta['beamModes'] if item['mode'] == stop), None)
        if not end_mode or not start_mode:
            self.status = Fill.STATUS.NORAMP
        else:

            # It's nice to have some context just before the ramp starts
            # Also, it's quite inexact
            # start_t = end_mode['startTime']
            start_t = start_mode['startTime']
            end_t = end_mode['endTime']
        
        # The actual data fetching
        if non_aligned_var:
            lg.log("fetching: {}...".format(", ".join(non_aligned_var)), end=" ")
            timber_vars = [self.timber_var_map[var] for var in non_aligned_var]
            data = self.db.get(timber_vars, start_t, end_t)
            for name, timber_var in self.timber_var_map.items():
                if not type(timber_var) == list and timber_var in data:
                    self.data[name] = Fill.Variable(data[timber_var])
            lg.log("done!", log_level=LogLevel.success, module_prestring=False)

        for var in aligned_var:
            lg.log("fetching aligned: {}...".format(var), end=' ')
            timber_vars = self.timber_var_map[var]
            data = self.db.getAligned(timber_vars, start_t, end_t)

            if len(data) == 0:
                raise Exception("No data found")
            xdata = data.pop('timestamps')
            ydata = np.stack((data[y] for y in data))
            self.data[var] = Fill.Variable((xdata, ydata))
            lg.log("done!", log_level=LogLevel.success, module_prestring=False)

        if cache:
            self.cache()

    def pack(self):
        return {
            'nbr' : self.nbr,
            'data' : self.data,
            'meta' : self.meta,
            'status' : self.status,
            'time_correction' : self.time_correction
        }

    def unpack(self, dump):
        self.nbr = dump['nbr']
        self.data = dump['data']
        self.meta = dump['meta']
        self.status = dump['status']
        self.time_correction = dump['time_correction'] if 'time_correction' in dump else self.time_correction

    def clear_cache(self):
        lg.log('clearing cache {}'.format(self.nbr))
        open(store_file_for_fill(self.nbr), 'wb')

    def cache(self):
        lg.log('caching {}'.format(self.nbr))
        with open(store_file_for_fill(self.nbr), 'wb') as f:
            pickle.dump(self.pack(), f)

    def load_cache(self):
        lg.log('loading {}'.format(self.nbr))
        with open(store_file_for_fill(self.nbr), 'rb') as f:
            self.unpack(pickle.loads(f.read()))

    # Data access methods
    def intensity(self):
        return self.data['intensity_b{}'.format(self.beam)]
    def blm_ir3(self):
        return self.data['synch_coll_b{}'.format(self.beam)]
    def blm_ir7(self):
        return self.data['beta_coll_b{}'.format(self.beam)]
    def motor_ir3(self):
        return self.data['motor_ir3_b{}'.format(self.beam)]
    def motor_ir7(self):
        return self.data['motor_ir7_b{}'.format(self.beam)]
    def motor_start(self):
        return self.data['motor_start']
    def abort_gap(self):
        return self.data['abort_gap_int_b{}'.format(self.beam)]
    def energy(self):
        return self.data['energy']

    def oml(self, aligned=False):
        if aligned:
            return self.data['A_synch_coll_b1']
        else:
            return self.data['synch_coll_b1']
    #### 


    ### Operations that affect the data somehow
    def normalize_intensity(self):
        var = 'intensity_b{}'.format(self.beam)
        self.data[var].y = self.data[var].y/np.max(self.data[var].y)

    def offset_time(self, align_mode="time_corr", t=0):
        # start, end = self.OML_period()
        # t = self.blm_ir3().x[start]

        if align_mode == "ramp":
            t = self.blm_ir3().x[0]
        elif align_mode == "ramp_mode":
            t = self.data['ramp_mode'].x[1]
        elif align_mode == "energy":
            t = self.energy().x[0]
        elif align_mode == "peak":
            # align w.r.t. beam 1 x peak
            beam = self.beam
            self.beam = 1
            # Previously we did alignment based of start of ramp which normally put the peak
            # at ~10 s. As we now align w.r.t. the peak, we move it 10 extra seconds so we
            # get similar results as before
            t = self.OML_peak()['t'] - 10 
            self.beam = beam
        elif align_mode == "time_corr":
            t = self.time_correction
        elif align_mode == "meas_time":
            ts = np.empty(self.motor_start().x.size)
            for i, v in enumerate(zip(*self.motor_start())):
                ts[i] = v[0] - v[1]*1e-9
            # t = stats.mode(ts)[0]
            t = np.median(ts)
            # print(t, self.motor_start().x[0])
        elif align_mode == "manual":
            t = self._timeshift = t
        else:
            raise Exception("align_mode does not exist")

        for v in self.data:
            self.data[v].x -= t
        self._timeshift = t

    def beta_coll_merge(self):
        # Is made obsolute 
        return

        var = 'beta_coll_b{}'.format(self.beam)
        if self.data[var].y.shape[0] == 3:
            self.data[var].y = np.sum(self.data[var].y, axis=0)

    ### Data queries
    def has_off_momentum_loss(self):
        """ Just some arbitrary test to see if we have a spike """
        max_loss = max(self.blm_ir3().y)
        mean_loss = np.mean(self.blm_ir3().y)
        if max_loss < 0 or mean_loss < 0:
            return False
        return max_loss/mean_loss > 10

    def OML_period(self): 
        """ The period is defined as
                start: t = 0
                end: t = cross over point
            returns (start, end) indices
        """ 
        i_co = self.crossover_point()['i']
        i_start = self.blm_ir3().index_for_time(0)
        return [i_start, i_co]

    def OML_peak(self):
        """ The timestamp and index for the maximum OML peak
            returns {'i' : ..., 't' : ...}
        """
        i_peak = imax(self.blm_ir3().y)[1]
        return {'i' : i_peak, 't' : self.blm_ir3().x[i_peak]}

    def crossover_point(self):
        """ Look for the point after OML spike when transversal losses starts 
            to dominate the momentum losses 
        """
        x = np.union1d(self.blm_ir3().x, self.blm_ir7().x)
        x = x[(x > self.blm_ir3().x.min())*(x < self.blm_ir3().x.max())]
        x = x[(x > self.blm_ir7().x.min())*(x < self.blm_ir7().x.max())]

        blm_ir3y = interpolate.interp1d(*self.blm_ir3())(x)
        blm_ir7y = interpolate.interp1d(*self.blm_ir7())(x)

        i = imax(blm_ir3y)[1]
        while i < len(x) - 1 and blm_ir3y[i] > blm_ir7y[i]: i += 1
        return {'t' : x[i], 'i' : self.blm_ir3().index_for_time(x[i])}

## class Fill
################


################
## Plotting directlyrelated to fills
def plot(fill):
    lg.log('plotting {}'.format(fill.nbr))
    fig = plt.figure()
    intensity_axis  = fig.add_subplot(211)
    energy_axis = intensity_axis.twinx()
    blm_axis = intensity_axis.twinx()


    intensity_axis.plot(*fill.intensity(), color='b', zorder=10, linestyle='-', linewidth=1)
    intensity_axis.set_ylim([0.95, 1.005])
    intensity_axis.set_ylabel("Frac. beam int.")

    energy_axis.plot(*fill.energy(), color='black', zorder=5)
    energy_axis.set_ylabel("Energy (GeV)")

    oml = fill.blm_ir3()
    start, end = fill.OML_period()
    blm_axis.axvspan(oml.x[start], oml.x[end], facecolor='b', alpha=0.2)

    fig.subplots_adjust(right=0.75)
    blm_axis.spines['right'].set_position(('axes', 1.15))
    blm_axis.set_frame_on(True)
    blm_axis.plot(*oml, color='r', linestyle='--', zorder=1, label='TCP IR3')
    blm_axis.set_yscale('log')
    blm_axis.set_ylabel("Losses (Gy/s)")
    blm_axis.plot(0,0, color='b', linestyle='-', label='Beam int.')
    blm_axis.plot(0,0, color='black', label='Ramp energy')
    blm_axis.legend(loc='lower right')

    plt.title("Fill {} (beam {})".format(fill.nbr, fill.beam))

    agap_ax = fig.add_subplot(212, sharex=intensity_axis)
    ag = moving_average(fill.abort_gap().y, 10)
    agap_ax.plot(fill.abort_gap().x, ag, color='g')
    agap_ax.set_ylim([0, 2e10])
    agap_ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "{0:.2f}".format(x/1e10)))
    agap_ax.set_ylabel("Abort gap int. ($10^{10}$ p.)")
    agap_ax.set_xlabel("Time (s)")

    plt.show()


def plot_blm(fill):
    lg.log('loss plot {}'.format(fill.nbr))
    fig = plt.figure()
    intensity_axis = fig.add_subplot(111)
    energy_axis = intensity_axis.twinx()
    blm_axis = intensity_axis.twinx()

    intensity_axis.plot(*fill.intensity(), color='b', zorder=10, linestyle='-', linewidth=1)
    intensity_axis.set_ylim([0.95, 1.005])
    intensity_axis.set_ylabel("Frac. beam int.")
    intensity_axis.set_xlabel("Time (s)")

    energy_axis.plot(*fill.energy(), color='black', zorder=5)
    energy_axis.set_ylabel("Energy (GeV)")

    oml = fill.blm_ir3()
    start, end = fill.OML_period()

    fig.subplots_adjust(right=0.8)
    blm_axis.spines['right'].set_position(('axes', 1.15))
    blm_axis.set_frame_on(True)
    blm_axis.plot(*oml, color='r', linestyle='--', zorder=2, label='TCP IR3')
    blm_axis.plot(*fill.blm_ir7(), color='g', linestyle='--', zorder=1, label='TCP-C IR7')
    blm_axis.axvspan(oml.x[start], oml.x[end], facecolor='b', alpha=0.2)
    blm_axis.set_yscale('log')
    blm_axis.set_ylabel("Losses (Gy/s)")

    blm_axis.plot(0,0, color='b', linestyle='-', label='Beam int.')
    blm_axis.plot(0,0, color='black', label='Ramp energy')
    blm_axis.legend(loc='lower right')


    blm_axis.set_xlim(fill.blm_ir3().x[[start, end]] + np.array([-5, +40]))

    plt.title("Fill {} (beam {})".format(fill.nbr, fill.beam))
    plt.show()

def plot_2_blm(f1, f2):
    lg.log("plot two")
    if f1.beam == 2 and f2.beam == 1:
        f1, f2 = (f2, f1)

    fig, blm_axis = plt.subplots()

    oml = f1.blm_ir3()
    start, end = f1.OML_period()

    blm_axis.spines['right'].set_position(('axes', 1.15))
    blm_axis.set_frame_on(True)
    
    blm_axis.plot(*oml, color='r', linestyle='-', zorder=2, label='1: TCP IR3')
    blm_axis.plot(*f1.blm_ir7(), color='g', linestyle='--', zorder=4, label='1: TCP-C IR7')

    blm_axis.plot(*f2.blm_ir3(), color='darkorange', linestyle='-', zorder=3, label='2: TCP IR3')
    blm_axis.plot(*f2.blm_ir7(), color='olivedrab', linestyle='--', zorder=5, label='2: TCP-C IR7')

    blm_axis.axvspan(oml.x[start], oml.x[end], facecolor='b', alpha=0.2)
    # blm_axis.axvline(x=46, color='black', zorder=0)

    blm_axis.set_yscale('log')
    blm_axis.set_ylabel("Losses (Gy/s)")
    blm_axis.set_xlabel("t (s)")
    blm_axis.legend(loc='upper right')
    blm_axis.set_xlim(f1.blm_ir3().x[[start, end]] + np.array([-5, +40]))

    f1_label = f1.nbr if f1.nbr > 0 else "AGGREGATE"
    f2_label = f2.nbr if f2.nbr > 0 else "AGGREGATE"
    if f1_label == f2_label:
        plt.title("Fill {}, beam 1 and 2".format(f1_label))
    else:
        plt.title("Fills {}b{}, {}b{}".format(f1_label, f1.beam, f2_label, f2.beam))
    plt.show()

def plot_motor(fill):
    lg.log("plotting motors")
    fig, motor_ax = plt.subplots()
    lax = motor_ax.twinx()

    motor_ax.plot(fill.motor_ir3().x, fill.motor_ir3().y[0], color='r')
    motor_ax.plot(fill.motor_ir3().x, fill.motor_ir3().y[1], color='r')
    motor_ax.set_ylabel("Motor (mm)")
    motor_ax.set_xlabel("Time (s)")

    lax.plot(*fill.blm_ir3())
    lax.set_ylabel("BLM")
    lax.set_yscale("log")

    plt.title("Motor plot: fill {}".format(fill.nbr))
    plt.show()

def plot_energy_ramp(fill):
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    plt.title("Energy ramp and derivative around start of ramp for fill 5427")

    ax1.plot(*fill.energy())
    ax1.set_ylabel("Energy GeV")

    ax2 = fig.add_subplot(212, sharex=ax1)
    d_energy = np.gradient(fill.energy().y)/np.gradient(fill.energy().x)
    ax2.plot(fill.energy().x, d_energy)
    ax2.set_ylabel("∆ Energy GeV")

    plt.show()



###############################
###### Some classification tools

def evaluate_off_momentum_losses_in_fills(fills, save_file):
    open(save_file, 'w').close() # erasing file

    for i in fills:
        lg.log("evaluating %s..." % str(i))
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
        lg.log("--\n")

def fills_from_file(file, status_string='*'):
    fills = []
    for line in open(file, 'r'):
        contents = line.rstrip('\n').split()
        fill_nbr = int(contents[0])
        status = contents[1]

        if status_string == "*" or status.upper() == status_string.upper():
            fills.append(fill_nbr)
    lg.log("Found {} fills".format(len(fills)))
    return fills

def plot_from(file, plot_at_the_time=10, status_string='*'):
    fills = fills_from_file(file, status_string)
    n = 0
    for fill_nbr in fills:
        n += 1
        lg.log("evaluating %s" % fill_nbr)
        fill = Fill(fill_nbr)
        plot(fill)
        lg.log("--\n")
        if n % plot_at_the_time == 0:
            inp = input("draw {} more plots? (press 'q' to quit) ".format(plot_at_the_time))
            if inp == 'q':
                break

def edit_event(e, axvspan, fill, fill_list):
    if e.key == "right":
        fill_list[fill.nbr][fill.beam][0] += 1
    if e.key == "left":
        fill_list[fill.nbr][fill.beam][0] -= 1
    oml = fill.blm_ir3()
    
    with open(Fill.OML_period_file, 'wb') as f:
        pickle.dump(fill_list, f)
        # json.dump(fill_list, f, indent=4, separators=(',', ': '), sort_keys=True)

    # I literally have no idea how this works, reference
    # http://stackoverflow.com/questions/35551903/update-location-of-axvspan-with-matplotlib
    arr = axvspan.get_xy()
    start = oml.x[fill_list[fill.nbr][fill.beam][0]]
    end = oml.x[fill_list[fill.nbr][fill.beam][1]]
    arr[:, 0] = [start, start, end, end, start]
    axvspan.set_xy(arr)
    plt.draw()
    

def edit_spike_for_fills(fills, status_string='OML'):
    try:
        with open(Fill.OML_period_file, 'rb') as f:
            fill_list = pickle.loads(f.read())
    except:
        fill_list = {}

    for nbr in fills:
        fill = Fill(nbr, fetch=False)
        fill.fetch()

        oml = fill.blm_ir3()
        start, end = fill.OML_period()
        try:
            fill_list[nbr][fill.beam] = [start, end]
        except:
            fill_list[nbr] = {fill.beam : [start, end]}

        fig, ax = plt.subplots()
        ax.plot(*oml, c='r')
        ax.set_yscale('log')
        ax.set_xlabel("t (s)")
        axvspan = ax.axvspan(oml.x[start], oml.x[end], facecolor='b', alpha=0.2)
        ax.set_xlim(oml.x[[start, start] + np.array([-10, 50])])
        fig.canvas.mpl_connect('key_press_event', lambda event : edit_event(event, axvspan, fill, fill_list))
        plt.title("Fill {} (beam {})".format(fill.nbr, fill.beam))
        plt.show()
    
    with open(Fill.OML_period_file, 'rb') as f:
        lg.log(pickle.loads(f.read()))


def intensity_and_OML_pruning(file_in, file_out):
    raise Exception("Has not been updated since refactoring")
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
        # smin, smax = find_spike(oml_data.y) 
        smin, smax = find_oml_spike(fill) 
        dur = oml_data.x[smax] - oml_data.x[smin]

        if dur < 100 or dur > 300:
            wrongly_categorised += 1
            continue

        with open(file_out, 'a') as f:
            f.write("{}\tOML\n".format(str(fill.nbr)))

    removed = int(round(float(low_intensity + wrongly_categorised)/len(fills) * 100))
    lg.log('discarded:')
    lg.log('\tintensity {}'.format(low_intensity))
    lg.log('\tmis-categorised {}'.format(wrongly_categorised))
    lg.log('\ttotal {}%'.format(removed))

def classify_time_alignment(fill_list):
    """ Classifying time offset with respect to TrimEditor defined motor variable
    """

    trimData = np.loadtxt("data/TCP.C6L7.B1--ramp.csv", delimiter=',', skiprows=2).transpose()
    trim = Fill.Variable((trimData[4], trimData[5]))
    
    for nbr in fill_list:
        fill = Fill(nbr, fetch=False)
        fill.fetch(forced=False)
        # fill.offset_time("peak")

        fd = Fill.Variable((np.array(fill.motor_ir7().x), np.array(fill.motor_ir7().y[1])))
        fd.x -= fd.x[0]

        # Fitting
        ip_trimx = interpolate.interp1d(trim.y, trim.x)(fd.y)
        func = lambda data, trim, c: np.sum((data + c - trim)**2)
        res = optimize.least_squares(func, 0.0, args=(fd.x, ip_trimx), jac='3-point')
        # print(res)

        fill.time_correction = fill.motor_ir7().x[0] + res.x[0]
        fill.cache()


def recache_fills(fill_list):
    for nbr in fill_list:
        for b in (1, 2):
            fill = Fill(nbr, beam=b, fetch=False)
            if b == 1: fill.clear_cache()
            fill.fetch(forced=(b == 1), cache=True)
        
