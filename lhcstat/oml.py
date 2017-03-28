import os
import pytimber
import matplotlib.pyplot as plt
import math
import ast
import pickle
import numpy as np
import bisect

from matplotlib.ticker import FuncFormatter

from scipy import stats
from scipy import interpolate


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

        def index_for_time(self, time):
            """ Helper function: find the index in the variable corresponding to the given timestamp """
            index = bisect.bisect_left(self.x, time)
            return index
    #### class Variable

    ## Variables in lists will be fetched aligned with each other
    BEAM1_VARIABLES = {
        'intensity_b1' : 'LHC.BCTFR.A6R4.B1:BEAM_INTENSITY',
        'beta_coll_b1' : ['BLMTI.06L7.B1E10_TCP.B6L7.B1:LOSS_RS09', 'BLMTI.06L7.B1E10_TCP.C6L7.B1:LOSS_RS09', 'BLMTI.06L7.B1E10_TCP.D6L7.B1:LOSS_RS09'],
        'synch_coll_b1' : 'BLMTI.06L3.B1I10_TCP.6L3.B1:LOSS_RS09',
        'abort_gap_int_b1' : 'LHC.BSRA.US45.B1:ABORT_GAP_INTENSITY',
    }

    BEAM2_VARIABLES = {
        'intensity_b2' : 'LHC.BCTFR.A6R4.B2:BEAM_INTENSITY',  
        'beta_coll_b2' : ['BLMTI.06R7.B2I10_TCP.B6R7.B2:LOSS_RS09', 'BLMTI.06R7.B2I10_TCP.C6R7.B2:LOSS_RS09', 'BLMTI.06R7.B2I10_TCP.D6R7.B2:LOSS_RS09'],
        'synch_coll_b2' : 'BLMTI.06R3.B2E10_TCP.6R3.B2:LOSS_RS09', 
        'abort_gap_int_b2' : 'LHC.BSRA.US45.B2:ABORT_GAP_INTENSITY'
    }

    COMB_VARIABLES = {
        'energy' : 'LHC.BOFSU:OFSU_ENERGY',
    }

    db = pytimber.LoggingDB()

    def __init__(self, nbr, fetch=True, beam=2):
        if not beam in (1, 2):
            raise Exception("Beam can only be 1 or 2")
        self.nbr = nbr
        self.beam = beam
        self.data = {}
        self.meta = {}
        self.status = Fill.STATUS.OK
        if self.beam == 1:
            self.timber_var_map = {**self.BEAM1_VARIABLES, **self.COMB_VARIABLES}
        else:
            self.timber_var_map = {**self.BEAM2_VARIABLES, **self.COMB_VARIABLES}

        if fetch:
            self.fetch()
            self.normalize_intensity()
            self.offset_time()


    def fetch(self, forced=False, cache=True):
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
        preramp = next((item for item in self.meta['beamModes'] if item['mode'] == 'PRERAMP'), None)
        ramp = next((item for item in self.meta['beamModes'] if item['mode'] == 'RAMP'), None)
        if not ramp or not preramp:
            self.status = Fill.STATUS.NORAMP
        else:
            start_t = ramp['startTime']
            end_t = ramp['endTime']
        
        # The actual data fetching
        if non_aligned_var:
            print("fetching: {}...".format(", ".join(non_aligned_var)), end=" ")
            timber_vars = [self.timber_var_map[var] for var in non_aligned_var]
            data = self.db.get(timber_vars, start_t, end_t)
            for name, timber_var in self.timber_var_map.items():
                if not type(timber_var) == list and timber_var in data:
                    self.data[name] = Fill.Variable(data[timber_var])
            print("done!")

        for var in aligned_var:
            print("fetching aligned: {}...".format(var), end=' ')
            timber_vars = self.timber_var_map[var]
            data = self.db.getAligned(timber_vars, start_t, end_t)

            xdata = data.pop('timestamps')
            ydata = np.stack((data[y] for y in data))
            self.data[var] = Fill.Variable((xdata, ydata))
            print("done!")

        if cache:
            self.cache()
        

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

    def cache(self):
        print('caching {}'.format(self.nbr))
        with open(store_file_for_fill(self.nbr), 'wb') as f:
            pickle.dump(self.pack(), f)

    def load_cache(self):
        print('loading {}'.format(self.nbr))
        with open(store_file_for_fill(self.nbr), 'rb') as f:
            self.unpack(pickle.loads(f.read()))

    # Data access methods
    def intensity(self):
        return self.data['intensity_b{}'.format(self.beam)]
    def blm_ir3(self):
        return self.data['synch_coll_b{}'.format(self.beam)]
    def blm_ir7(self):
        return self.data['beta_coll_b{}'.format(self.beam)]
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

    def offset_time(self, by="oml"):
        if by == "oml":
            start, end = self.OML_period()
            t = self.blm_ir3().x[start]
        elif by == "ramp":
            raise Exception("deprecated")
            istart = find_start_of_ramp(self.data["energy"])
            t = self.data["energy"].x[istart]
        else:
            raise Exception("don't recognize offset type")

        for v in self.data:
            self.data[v].x -= t

    def beta_coll_merge(self):
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
                start: when OML first is noticed
                end: when the losses (and the derivative) has gone below a certain threshold
            returns (start, end) indices
        """ 
        data = self.blm_ir3().y
        vpeak, ipeak = imax(data)
        start = end = ipeak
    
        ddata = np.abs(np.gradient(data))
        threshold = 1e-7
        # with open('fills/spikes.dat', 'rb') as f:
            # spike_list = pickle.loads(f.read())
            # if self.nbr in spike_list:
                # return spike_list[self.nbr]
        
        found_start = found_end = False
        while not found_start and start > 0:
            start -= 1
            # if ddata[start] < threshold and data[start] < vpeak/5e1:
            if ddata[start] < 1e-5:
                found_start = True
        
        while not found_end and end < len(ddata) - 1:
            end += 1
            if ddata[end] < threshold and data[end] < vpeak/1e2:
                found_end = True
        
        return (start, end)

    def crossover_point(self):
        """ Look for the point after OML spike when transversal losses starts 
            to dominate the momentum losses 
    
            Note: this should be used with aligned data """
        self.beta_coll_merge()
        oml = self.blm_ir3() # off momentum loss
        tl = self.blm_ir7() # transversal losses
        vpeak, ipeak = imax(oml.y)
    
        i = ipeak
        while oml.y[i] > tl.y[i]:
            i += 1
        return {'i' : i, 't' : oml.x[i]}

## class Fill
################


################
## Plotting directlyrelated to fills
def plot(fill):
    print('plotting {}'.format(fill.nbr))
    fig = plt.figure()
    intensity_axis  = fig.add_subplot(211)
    energy_axis = intensity_axis.twinx()
    blm_axis = intensity_axis.twinx()


    intensity_axis.plot(*fill.intensity(), color='b', zorder=10, linestyle='-', linewidth=1)
    intensity_axis.set_ylim([0.95, 1.005])
    intensity_axis.set_ylabel("Beam Intensity")

    energy_axis.plot(*fill.energy(), color='black', zorder=5)
    energy_axis.set_ylabel("Energy")

    oml = fill.blm_ir3()
    start, end = fill.OML_period()
    blm_axis.axvspan(oml.x[start], oml.x[end], facecolor='b', alpha=0.2)

    fig.subplots_adjust(right=0.8)
    blm_axis.spines['right'].set_position(('axes', 1.15))
    blm_axis.set_frame_on(True)
    blm_axis.plot(*oml, color='r', linestyle='--', zorder=1)
    blm_axis.set_yscale('log')
    blm_axis.set_ylabel("Losses")


    plt.title("Fill {}".format(fill.nbr))

    agap_ax = fig.add_subplot(212, sharex=intensity_axis)
    agap_ax.plot(fill.abort_gap().x, moving_average(fill.abort_gap().y, 10), color='g')
    agap_ax.set_ylabel("Abort gap intensity")
    agap_ax.set_xlabel("Time (s)")

    plt.show()


def plot_blm(fill):
    print('loss plot {}'.format(fill.nbr))
    fig = plt.figure()
    intensity_axis = fig.add_subplot(111)
    energy_axis = intensity_axis.twinx()
    blm_axis = intensity_axis.twinx()

    intensity_axis.plot(*fill.intensity(), color='b', zorder=10, linestyle='-', linewidth=1)
    intensity_axis.set_ylim([0.95, 1.005])
    intensity_axis.set_ylabel("Beam Intensity")
    intensity_axis.set_xlabel("Time (s)")

    energy_axis.plot(*fill.energy(), color='black', zorder=5)
    energy_axis.set_ylabel("Energy")

    oml = fill.blm_ir3()
    start, end = fill.OML_period()

    fig.subplots_adjust(right=0.8)
    blm_axis.spines['right'].set_position(('axes', 1.15))
    blm_axis.set_frame_on(True)
    blm_axis.plot(*oml, color='r', linestyle='--', zorder=2, label='momentum loss')
    blm_axis.plot(*fill.blm_ir7(), color='g', linestyle='--', zorder=1, label='transversal loss')
    blm_axis.axvspan(oml.x[start], oml.x[end], facecolor='b', alpha=0.2)
    blm_axis.set_yscale('log')
    blm_axis.set_ylabel("Losses")
    blm_axis.legend(loc='lower right')

    plt.title("Fill {}".format(fill.nbr))
    plt.show()

def plot_motor(fill):
    raise Exception("need to fetch motor values -- not done right now")
    print("plotting motors")
    fig, motor_ax = plt.subplots()
    e_ax = motor_ax.twinx()

    motor_ax.plot(*fill.data['A_tcp_motor_ru'], color='r')
    motor_ax.plot(*fill.data['A_tcp_motor_lu'], color='r')
    motor_ax.set_ylabel("Motor (mm σ)")
    motor_ax.set_xlabel("Time (s)")

    e_ax.plot(*fill.data['energy'])
    e_ax.set_ylabel("Energy (MeV)")

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
        #     fill.crop_ramp()
        plot(fill)
        print("--\n")
        if n % plot_at_the_time == 0:
            inp = input("draw {} more plots? (press 'q' to quit) ".format(plot_at_the_time))
            if inp == 'q':
                break

def edit_event(e, axvspan, fill, fill_list):
    raise Exception("Has not been updated since refactoring")
    if e.key == "right":
        fill_list[fill.nbr][0] += 1
    if e.key == "left":
        fill_list[fill.nbr][0] -= 1
    oml = fill.data['synch_coll_b1']

    with open("fills/spikes.dat", 'wb') as f:
        pickle.dump(fill_list, f)

    # I literally have no idea how this works, reference
    # http://stackoverflow.com/questions/35551903/update-location-of-axvspan-with-matplotlib
    arr = axvspan.get_xy()
    start = oml.x[fill_list[fill.nbr][0]]
    end = oml.x[fill_list[fill.nbr][1]]
    arr[:, 0] = [start, start, end, end, start]
    axvspan.set_xy(arr)
    plt.draw()
    

def edit_spike_for_fills(file, status_string='OML'):
    raise Exception("Has not been updated since refactoring")
    fills = fills_from_file(file, status_string)
    fill_list = {}

    c = 0
    for nbr in fills:
        fill = Fill(nbr)

        oml = fill.data['synch_coll_b1']
        start, end = find_oml_spike(fill)
        fill_list[nbr] = [start, end]

        fig, ax = plt.subplots()
        ax.plot(*oml, c='r')
        ax.set_yscale('log')
        ax.set_xlabel("t (s)")
        axvspan = ax.axvspan(oml.x[start], oml.x[end], facecolor='b', alpha=0.2)
        ax.set_xlim([oml.x[start] - 20, oml.x[start] + 20])
        fig.canvas.mpl_connect('key_press_event', lambda event : edit_event(event, axvspan, fill, fill_list))
        plt.show()
    
    with open("fills/spikes.dat", 'rb') as f:
        print(pickle.loads(f.read()))


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
    print('discarded:')
    print('\tintensity {}'.format(low_intensity))
    print('\tmis-categorised {}'.format(wrongly_categorised))
    print('\ttotal {}%'.format(removed))

