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

    # For beam 2
    #   L  -> R
    #   B1 -> B2
    b1_variables = {
            }
    variables = {
        'intensity_b1' : 'LHC.BCTFR.A6R4.B1:BEAM_INTENSITY',
        # 'intensity_b2' : 'LHC.BCTFR.A6R4.B2:BEAM_INTENSITY',
        'energy' : 'LHC.BOFSU:OFSU_ENERGY',
        'beta_coll_b_b1' : 'BLMTI.06L7.B1E10_TCP.B6L7.B1:LOSS_RS09',
        'beta_coll_c_b1' : 'BLMTI.06L7.B1E10_TCP.C6L7.B1:LOSS_RS09',
        'beta_coll_d_b1' : 'BLMTI.06L7.B1E10_TCP.D6L7.B1:LOSS_RS09',
        'synch_coll_b1' : 'BLMTI.06L3.B1I10_TCP.6L3.B1:LOSS_RS09',
        'abort_gap_int_b1' : 'LHC.BSRA.US45.B1:ABORT_GAP_INTENSITY',

        # These might be for B2
        # 'BLMTI.06R3.B2E10_TCP.6R3.B2:LOSS_RS09',
        # 'BLMTI.06R7.B2I10_TCP.B6R7.B2:LOSS_RS09',
        # 'BLMTI.06R7.B2I10_TCP.C6R7.B2:LOSS_RS09',
        # 'BLMTI.06R7.B2I10_TCP.D6R7.B2:LOSS_RS09',
    }

    aligned_variables = {
        ## ALIGNED DATA
        ## Needs to have a A_ prefix
        'A_beta_coll_b_b1' : 'BLMTI.06L7.B1E10_TCP.B6L7.B1:LOSS_RS09',
        'A_beta_coll_c_b1' : 'BLMTI.06L7.B1E10_TCP.C6L7.B1:LOSS_RS09',
        'A_beta_coll_d_b1' : 'BLMTI.06L7.B1E10_TCP.D6L7.B1:LOSS_RS09',
        'A_synch_coll_b1' : 'BLMTI.06L3.B1I10_TCP.6L3.B1:LOSS_RS09',
        #
        # Motors
        'A_tcp_motor_ld' : 'TCP.6L3.B1:MEAS_MOTOR_LD',
        'A_tcp_motor_lu' : 'TCP.6L3.B1:MEAS_MOTOR_LU',
         'A_tcp_motor_rd' : 'TCP.6L3.B1:MEAS_MOTOR_RD',
         'A_tcp_motor_ru' : 'TCP.6L3.B1:MEAS_MOTOR_RU',
    }
    all_variables = {**variables, **aligned_variables}

    db = pytimber.LoggingDB()

    def __init__(self, nbr, fetch=True, beam=1):
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
        self.fetch_2(forced, cache)
        return

        to_fetch = self.variables
        fetch_aligned = forced
        if not forced and os.path.isfile(store_file_for_fill(self.nbr)):
            self.load_cache()

            non_cached_variables = []
            for v in self.all_variables:
                if not v in self.data.keys():
                    if v.startswith("A_"): # We have some unfetched aligned variable -> we need to refetch all aligned
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
            print("\tfetch", to_fetch)
            data = self.db.get([self.all_variables[v] for v in to_fetch], start_t, end_t)
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
            self.cache()

    def fetch_2(self, forced=False, cache=True):
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

##### class Fill
################


## PLOTTING OF FILL DATA
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
    raise Exception("Has not been updated")
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
    raise Exception("Has not been updated")
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


def find_start_of_ramp(energy): 
    raise Exception("Deprecated")
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


def find_oml_spike(fill): 
    raise Exception("Deprecated")
    data = fill.data['synch_coll_b1'].y
    vpeak, ipeak = imax(data)
    start = end = ipeak

    ddata = np.abs(np.gradient(data))
    threshold = 1e-7
    with open('fills/spikes.dat', 'rb') as f:
        spike_list = pickle.loads(f.read())
        if fill.nbr in spike_list:
            return spike_list[fill.nbr]
    
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


def find_crossover_point(fill):
    """ Look for the point after OML spike when transversal losses starts 
        to dominate the momentum losses 

        Note: this should be used with aligned data """
    raise Exception("Deprecated")
    fill.beta_coll_merge()
    oml = fill.data["synch_coll_b1"] # off-momentum-losses
    tl = fill.data["A_beta_coll_b1"]   # transversal-losses
    vpeak, ipeak = imax(oml.y)

    i = ipeak
    while oml.y[i] > tl.y[i]:
        i += 1
    return {'i' : i, 't' : oml.x[i]}


def merge_variable(fills, var):
    raise Exception("moved")
    xmax, xmin = 0, 0
    for fill in fills:
        xmin = min(min(fill.data[var].x), xmin)
        xmax = max(max(fill.data[var].x), xmax)

    # The x-series
    newx = np.arange(xmin, xmax + 1, 1, dtype=np.int) # using int's to maybe make comparisons more safe

    # Will store the average values of y
    avgy = np.zeros(newx.size, dtype=np.float)
    norm = np.zeros(newx.size, dtype=np.float) # we use this to basically take the average value at each timestamp

    # Will store maximum and minimum values of y
    maxy = np.zeros(newx.size, dtype=np.float)
    miny = 1.0e10*np.ones(newx.size, dtype=np.float)

    for fill in fills:
        x_i = fill.data[var].x.astype(np.int)
        index_map = np.searchsorted(newx, x_i)

        np.add.at(avgy, index_map, fill.data[var].y)
        np.add.at(norm, index_map, 1.0)

        np.maximum.at(maxy, index_map, fill.data[var].y)
        np.minimum.at(miny, index_map, fill.data[var].y)

    # delete empty values (where no data was written)
    to_delete = np.where(norm < 0.5)
    newx = np.delete(newx, to_delete)
    maxy = np.delete(maxy, to_delete)
    miny = np.delete(miny, to_delete)
    norm = np.delete(norm, to_delete)
    avgy = np.delete(avgy, to_delete)
    avgy /= norm

    newx = newx.astype(dtype=float)

    return {
            "average" : Fill.Variable([newx, avgy]),
            "max" : Fill.Variable([newx, maxy]),
            "min" : Fill.Variable([newx, miny])
            }


def aggregate_fill(fill_list=[], from_cache=False):
    """ Create an aggregate fill of the fill ids. If reading from cache, 
        the fill_list can be empty. """

    raise Exception("moved")
    if not fill_list and not from_cache:
        raise Exception("'fill_list' can't be empty if not to read from cache'")

    aggr_id = 0
    fill = Fill(aggr_id, fetch=False)
    if from_cache:
        fill.load_cache()
        return fill

    fills = []
    for nbr in fill_list:
        fill = Fill(nbr)
        fill.beta_coll_merge()
        fills.append(fill)
        
    var = ["synch_coll_b1", "A_beta_coll_b1", "energy", "intensity_b1"]
    for v in var:
        merged = merge_variable(fills, v)
        fill.data[v] = merged["average"]
    return fill


def plot_aggregate_fill(fill_list):
    raise Exception("moved")

    # We can't reuse the 'aggregate_fill' function above, as we want to plot more
    # than for a normal fill (min, max, average)
    fills = []
    for nbr in fill_list:
        fill = Fill(nbr)
        fill.beta_coll_merge()
        fills.append(fill)
        
    var = ["synch_coll_b1", "A_beta_coll_b1", "energy"]
    merged = {}
    for v in var:
        merged[v] = merge_variable(fills, v)

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    e_ax = ax.twinx()

    v = "synch_coll_b1"
    ax.plot(*merged[v]["average"], zorder=10, color='red', label='OML')
    ax.fill_between(*merged[v]["max"], merged[v]["min"].y, color='red', alpha='0.3')

    v = "A_beta_coll_b1"
    ax.plot(*merged[v]["average"], zorder=9, color='green', label='TL')
    ax.fill_between(*merged[v]["max"], merged[v]["min"].y, color='green', alpha='0.2')

    e_ax.plot(*merged["energy"]["average"], zorder=5, color='black', label='energy')
    e_ax.set_ylabel("Reference energy (GeV)")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("BLM signal")
    ax.legend(loc="upper right")
    ax.set_xlim([-20, 60])
    # ax.axvspan(0.0, 13.55, facecolor='b', zorder=0, alpha=0.1)
    plt.title("Aggregate fill 2016")
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

def export_energy_ramp_to_sixtrack(file_out='ramp.txt', fill=5433, interpolation_type='cubic'):
    raise Exception("Deprecated")
    f = Fill(fill)

    file_out = "calc/{}_{}".format(fill, file_out)
    print("Ramp export: {} interpolation".format(interpolation_type))
    print("Writing to file '{}'".format(file_out))

    freq = 11245 # laps / s
    with open(file_out, 'w') as file:
        x = f.data['energy'].x[0:300] # at least 300 s, should be enough
        y = f.data['energy'].y[0:300]
        rescaled_x = [freq*(i - x[0]) for i in x]
        delta = int(round(rescaled_x[-1] - rescaled_x[0]))

        new_x = np.linspace(0, delta, delta)
        new_y = interpolate.interp1d(rescaled_x, y, kind=interpolation_type)(new_x)
        for i, e in enumerate(new_y):
            scaled_e = e*1e3 # Sixtrack uses MeV
            turnnbr = int(round(new_x[i])) + 1 # 1 indexed
            file.write('{} {}\n'.format(turnnbr, scaled_e))

def export_momentum_tcp_squeeze(file_out='motor_tcp.txt', fill=5433, interpolation_type='cubic'):
    raise Exception("Deprecated")
    f = Fill(fill)

    file_out = "calc/{}_{}".format(fill, file_out)
    print("TCP movement: {} interpolation".format(interpolation_type))
    print("Writing to file '{}'".format(file_out))

    ispike = find_start_of_ramp(f.data['energy'])
    tspike = f.data['energy'].x[ispike]
    istart = bisect.bisect_left(f.data['A_tcp_motor_ld'].x, tspike)
    iend = istart + 300 # at least 300 s

    dispersion = -2.07e3 # in mm
    freq = 11245 # laps / s
    with open(file_out, 'w') as file:
        # the data is aligned, so we can use the same x for both y
        x = f.data['A_tcp_motor_ld'].x[istart:iend]
        y_low = f.data['A_tcp_motor_ld'].y[istart:iend]
        y_high = f.data['A_tcp_motor_rd'].y[istart:iend]
        rescaled_x = [freq*(i - x[0]) for i in x]
        delta = int(round(rescaled_x[-1] - rescaled_x[0]))

        new_x = np.linspace(0, delta, delta)
        new_y_low = interpolate.interp1d(rescaled_x, y_low, kind=interpolation_type)(new_x)
        new_y_high = interpolate.interp1d(rescaled_x, y_high, kind=interpolation_type)(new_x)
        for i, m_low in enumerate(new_y_low):
            m_high = new_y_high[i]
            turnnbr = int(round(new_x[i])) + 1 # 1 indexed
            file.write('{} {} {}\n'.format(turnnbr, m_low/dispersion, m_high/dispersion))

#################
## Statistics

def draw_histogram(title, data, binsize, xlabel='', ylabel='', color='b'):
    raise Exception("moved")
    maxbin = max(data) + binsize
    minbin = min(data)
    bins = np.arange(minbin, maxbin, binsize)
    fig, ax = plt.subplots()
    ax.hist(data, bins=bins, color=color, edgecolor='black')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()


def intensity_histogram(file):
    raise Exception("moved")
    fills = fills_from_file(file, "OML")
    intensities = []
    for nbr in fills:
        fill = Fill(nbr, fetch=False)
        fill.fetch()
        # fill.crop_ramp()

        intensities.append(max(fill.data['intensity_b1'].y))

    draw_histogram('Intensity for {}'.format(file), intensities, 1e13, "Intensity", "Count")

def spike_duration_histogram(file):
    raise Exception("moved")
    fills = fills_from_file(file, "OML")
    outliers = []
    durations = []
    for nbr in fills:
        fill = Fill(nbr)

        # losses = fill.data['synch_coll_b1'].y
        # spike_start, spike_end = find_spike(np.array(losses))
        spike_start, spike_end = find_oml_spike(fill)
        d = fill.data['synch_coll_b1'].x[spike_end] - fill.data['synch_coll_b1'].x[spike_start]
        if d < 70 or d > 300:
            outliers.append(nbr)
        durations.append(d)

    draw_histogram('Spike duration for {}'.format(file), durations, 10, 'Seconds', 'Count')
    return outliers

def oml_dom_duration_histogram(file):
    raise Exception("moved")
    """ Histogram on what durations OML dominate the transversal losses """
    fills = fills_from_file(file, "OML")
    durations = []
    outliers = []
    for nbr in fills:
        fill = Fill(nbr)
        
        crossover = find_crossover_point(fill)
        if crossover['t'] > 40:
            outliers.append(nbr)
        else:
            durations.append(crossover["t"])

    draw_histogram("Duration OML > transversal losses from '{}'".format(file), durations, 1, 'Duration (s) after spike with OML > TM', 'Count')
    return outliers

def fills_bar_graph(file):
    raise Exception("moved")
    """ Bar graph displaying the 'time till max spike' + 'time where OML > TM' for all fills in file """
    fills = fills_from_file(file, "OML")
    oml_dom_duration = []
    time_till_spike = []
    total_duration = []
    for nbr in fills:
        fill = Fill(nbr)

        crossover = find_crossover_point(fill)
        if crossover['t'] < 40: # why am I using this?
            oml = fill.data['synch_coll_b1']
            ispike = np.argmax(oml.y)
            tspike = oml.x[ispike]
            time_till_spike.append(tspike)

            oml_dom_duration.append(crossover["t"] - tspike)

            total_duration.append(crossover['t'])

    fig, ax = plt.subplots()
    ax.bar(range(len(oml_dom_duration)), time_till_spike, label='Time (s) till spike')
    ax.bar(range(len(oml_dom_duration)), oml_dom_duration, bottom=time_till_spike, label='Time (s) where OML > TM')
    ax.legend(loc="upper right")
    ax.set_xlabel("Range id")
    ax.set_ylabel("Duration (s)")
    plt.show()


def max_spike_histogram(file):
    raise Exception("moved")
    fills = fills_from_file(file, "OML")
    spike_time = []

    for nbr in fills:
        fill = Fill(nbr)

        losses = fill.data['synch_coll_b1']
        ispike = np.argmax(losses.y)
        spike_time.append(losses.x[ispike])

    draw_histogram('Max spike event for {}'.format(file), spike_time, 1.0, 'Delta t (s) from start of ramp till spike', 'Count')

def beta_vs_synch_blm(file):
    raise Exception("moved")
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

        # smin, smax = find_spike(fill.data['synch_coll_b1'].y) 
        smin, smax = find_oml_spike(fill) 
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
    raise Exception("moved")
    fills = fills_from_file(file, "OML")
    intensity = []
    mean_loss = []
    max_loss = []
    discarded = 0
    for nbr in fills:
        fill = Fill(nbr, False)
        fill.fetch()
        # smin, smax = find_spike(fill.data['synch_coll_b1'].y) 
        smin, smax = find_oml_spike(fill) 
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
    raise Exception("moved")
    fills = fills_from_file(file, "OML")
    abort_gap = []
    average_loss = []
    max_loss = []
    for nbr in fills:
        fill = Fill(nbr, False)
        fill.fetch()
        smin, smax = find_oml_spike(fill) 

        # Only looking until t_co instead -- will not affect max
        smax = find_crossover_point(fill)['i']

        tmax = fill.data['synch_coll_b1'].x[smax]
        tmin = fill.data['synch_coll_b1'].x[smin]

        # tmax = find_crossover_point(fill)['t']

        ag_average = moving_average(fill.data['abort_gap_int_b1'].y, 5)
        agmin = fill.data['abort_gap_int_b1'].index_for_time(tmin)
        agmax = fill.data['abort_gap_int_b1'].index_for_time(tmax)

        ssubset = fill.data['synch_coll_b1'].y[smin:smax]

        average_loss.append(np.average(ssubset))
        max_loss.append(max(ssubset))
        abort_gap.append(ag_average[agmin] - ag_average[agmax])


    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, sharey=ax1) 

    # fig1, ax1 = plt.subplots()
    ax1.set_xlabel("Average BLM")
    ax1.set_ylabel("∆ abort gap intensity")
    ax1.scatter(average_loss, abort_gap, color='b', label='average')
    ax1.set_xlim([0, 1.1*max(average_loss)])
    ax1.set_ylim([0, 1.1*max(abort_gap)])

    xval = [0, 1]
    slope, intercept, r_value, p_value, std_err = stats.linregress(average_loss, abort_gap)
    print("Average fit")
    print("\tk  ={:>10.3E}\n\tm  ={:>10.3E}\n\tr  ={:>10.7f}\n\tp  ={:>10.3E}\n\te^2={:>10.3E}".format(slope, intercept, r_value, p_value, std_err))
    yfit = [slope*x + intercept for x in xval]
    ax1.plot(xval, yfit, color='gray')

    ax1.legend(loc="lower right")

    # fig2, ax2 = plt.subplots()
    ax2.set_xlabel("Max BLM")
    ax2.scatter(max_loss, abort_gap, color='r', label='max')
    ax2.set_xlim([0, 1.1*max(max_loss)])
    ax2.legend(loc="lower right")

    slope, intercept, r_value, p_value, std_err = stats.linregress(max_loss, abort_gap)
    print("Max fit")
    print("\tk  ={:>10.3E}\n\tm  ={:>10.3E}\n\tr  ={:>10.7f}\n\tp  ={:>10.3E}\n\te^2={:>10.3E}".format(slope, intercept, r_value, p_value, std_err))
    yfit = [slope*x + intercept for x in xval]
    ax2.plot(xval, yfit, color='gray')

    fig.suptitle("Correlation between abort gap intensity and BLM signal for TCP in IR3")
    plt.show()
    
def abort_gap_vs_BLM(file):
    raise Exception("moved")
    fills = fills_from_file(file, "OML")
    abort_gap = []
    smean_loss = []
    bmean_loss = []
    for nbr in fills:
        fill = Fill(nbr, False)
        fill.fetch()
        fill.beta_coll_merge()

        # smin, smax = find_spike(fill.data['synch_coll_b1'].y) 
        smin, smax = find_oml_spike(fill) 
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
#     fills = fills_from_file(file, "OML")
#     integrated_losses = []
#     for nbr in fills:
#         fill = Fill(nbr)

#         losses = np.array(fill.data['synch_coll_b1'][1])
#         spike, tail = find_spike(losses)
#         int_loss = 0.0
#         for i in range(spike, tail):
#             int_loss += losses[i]
#         integrated_losses.append(int_loss)

#     draw_histogram('Integrated losses', integrated_losses, 0.01)
#     return integrated_losses


# def max_loss_histogram(file):
#     fills = fills_from_file(file, "OML")
#     max_loss = []
#     for nbr in fills:
#         fill = Fill(nbr, fetch=False)
#         fill.fetch()
#         # fill.crop_ramp() # can throw
#         max_loss.append(max(fill.data['synch_coll_b1'][1]))

#     draw_histogram('Max losses', max_loss, 0.005)
#     return max_loss

# def find_spike(data):
#     raise Exception("should not be used")
#     # ddata = np.gradient(data)
#     spike_index = max(range(len(data)), key=data.__getitem__) # - 3 # just so we get more of the start of the spike
#     # spike_val = data[spike_index]

#     # Used 21 previously --> feels unreasonable
#     mov_average = moving_average(data, 5)
#     # dmov_average = np.gradient(mov_average)

#     # Use the moving average to calculate the tail
#     average_max = max(mov_average)
#     tail_threshold = average_max/10.0
#     tail_index = -1
#     for i, d in enumerate(mov_average[spike_index:]):
#         if d < tail_threshold:
#             tail_index = spike_index + i
#             break
#     else:
#         print("Spike: ", spike_index, tail_index)
#         raise Exception("did not find tail")

#     # Note that
#     #    spike_index < tail_index
#     return [spike_index, tail_index]

