import numpy as np
import pytimber
import os
import pickle
from collections import namedtuple
from scipy import interpolate

class CollTrimVMap:

    fillVars = np.array([
        "t",
        "MEAS_MOTOR_LD",
        "MEAS_MOTOR_RD",
        "MEAS_MOTOR_LU",
        "MEAS_MOTOR_RU",
        "MEAS_LIMIT_DUMP_INNER_LD",
        "MEAS_LIMIT_DUMP_INNER_LU",
        "MEAS_LIMIT_DUMP_INNER_RD",
        "MEAS_LIMIT_DUMP_INNER_RU",
        "MEAS_LIMIT_DUMP_OUTER_LD",
        "MEAS_LIMIT_DUMP_OUTER_LU",
        "MEAS_LIMIT_DUMP_OUTER_RD",
        "MEAS_LIMIT_DUMP_OUTER_RU",
        "MEAS_LIMIT_DUMP_INNER_GU",
        "MEAS_LIMIT_DUMP_INNER_GD",
        "MEAS_LIMIT_DUMP_OUTER_GU",
        "MEAS_LIMIT_DUMP_OUTER_GD",
        ])

    trimVars = np.array([
        "t",
        "left_downstream",
        "right_downstream",
        "left_upstream",
        "right_upstream",
        "dump_inner_left_downstream",
        "dump_inner_left_upstream",
        "dump_inner_right_downstream",
        "dump_inner_right_upstream",
        "dump_outer_left_downstream",
        "dump_outer_left_upstream",
        "dump_outer_right_downstream",
        "dump_outer_right_upstream",
        "dump_inner_gap_upstream",
        "dump_inner_gap_downstream",
        "dump_outer_gap_upstream",
        "dump_outer_gap_downstream",
        ])

    @classmethod
    def trimVarsForFill(clss, fill_vars):
        return clss.trimVars[np.in1d(clss.fillVars, fill_vars)]

    @classmethod
    def fillVarsForTrim(clss, trim_vars):
        return clss.fillVars[np.in1d(clss.trimVars, trim_vars)]

class CollFunction:
    def __init__(self, cname):
        self.cname = cname
        self.dmap = np.empty(0)
        self.data = np.empty(0)

    def get_series(self, names):
        """ Get data for given names. 
        
            E.g.  
                get_series(names=['t', 'left_upstream', 'left_downstream'])
            returns np.array with the respective columns.
        """
        m = np.empty(len(names), dtype=int)
        for i, n in enumerate(names):
            match = np.where(self.dmap == n)[0]
            if match.size == 0: raise Exception("no '{}' in dmap".format(n))
            m[i] = match[0]
        return self.data[m, :]
    
    def recalculate_time(self, func):
        """ Recalculate the time variable so that the 
            two functions can be compared per index basis.

            Does interpolation. Values outside the previous
            range will be set to 0.
        """
        t = self.get_series("t")[0]
        other_t = func.get_series("t")[0]
        new_data = np.empty((self.data.shape[0], other_t.size))
        for i, d in enumerate(self.data):
            if self.dmap[i] == "t":
                new_data[i] = other_t
            new_data[i] = interpolate.interp1d(t, d, bounds_error=False, fill_value="extrapolate")(other_t)
        self.data = new_data
        

class CollTrimFunction(CollFunction):
    def __init__(self, csv_header, csv_data_block):
        """ Parse the data as received from the csv file. 
            
            We assume:
                1. an even amount of columns
                2. the sample time _is the same_ for all motors
                3. the data is distributed as 't, motor1, t, motor2, ...'
        """

        if not csv_header.shape[0] == csv_data_block.shape[1] or not len(csv_header) % 2 == 0:
            raise Exception("illegal data: unexpected shape")
        cnames = np.unique([n.split('/')[0] for n in csv_header])
        if cnames.size > 1:
            raise Exception("more than one collimator in header")

        super(CollTrimFunction, self).__init__(cnames[0])

        # Use the first column as the global time, remove the others
        mtypes = np.array([n.split('#')[-1] for n in csv_header])
        mtypes[0] = 't'
        delete = np.arange(2, mtypes.size, 2)
        self.dmap = np.delete(mtypes, delete)
        self.data = np.delete(csv_data_block, delete, axis=1).transpose()

class CollFillFunction(CollFunction):
    db = pytimber.LoggingDB()
    cache_dir = "cache"

    def __init__(self, fill_nbr, cname):
        super(CollFillFunction, self).__init__(cname)
        self.fill_nbr = fill_nbr

    def fetch(self, forced=False, cache=True):
        """ Fetch all logged motor functions for collimator 'cname'
            and fill 'nbr'.
        """
        print("loading collimator {}... ".format(self.cname), end=" ")
        if not forced and os.path.isfile(self.cache_file()):
            print("from cache")
            with open(self.cache_file(), 'rb') as f:
                self.unpack(pickle.loads(f.read()))
            return

        print("fetching", end=" ")
        # Fill data
        meta = self.db.getLHCFillData(self.fill_nbr)
        ramp_mode = next((item for item in meta['beamModes'] if item['mode'] == "RAMP"), None)
        if not ramp_mode:
            raise Exception("did not find ramp beam mode in fill")
        start = ramp_mode['startTime']
        end = ramp_mode['endTime']
        
        # Collimator motor data
        variables = ["{}:{}".format(self.cname, v) for v in CollTrimVMap.fillVars[CollTrimVMap.fillVars != "t"]]
        variables.append("{}:MEAS_PROFILE_TIME".format(self.cname))
        print(variables)
        response = self.db.getAligned(variables, start, end)
        dmap = []
        data = []
        print(response.keys())
        for key in response:
            if "PROFILE" in key: 
                i_align = np.argmin(response[key])
                continue
            elif key == "timestamps":
                dmap.append('t')
            else:
                dmap.append(key.split(":")[1])
            data.append(response[key])
        self.dmap = np.array(dmap)
        self.data = np.array(data)

        # align w.r.t. profile function
        t = np.where(self.dmap == 't')[0]
        self.data[t] -= self.data[t, i_align]

        if cache:
            with open(self.cache_file(), 'wb') as f:
                pickle.dump(self.pack(), f)

    def pack(self):
        return {
            'dmap' : self.dmap,
            'data' : self.data
            }
    
    def unpack(self, dump):
        self.dmap = dump['dmap']
        self.data = dump['data']

    def cache_file(self):
        return "{}/{}_{}.dat".format(self.cache_dir, self.fill_nbr, self.cname)
        
    
def cfunc_from_trimfile(filepath):
    colls = []
    with open(filepath, 'r') as f:
        header = np.array(f.readline().strip().split(','))
        f.readline()
        data = np.loadtxt(f, delimiter=',')

        m = np.array([not "warning" in h for h in header], dtype=bool)
        header = header[m]
        data = data[:, m]
    
        i = 0
        while i < len(header):
            cname = header[i].split('/')[0]

            start = i
            while i < len(header):
                cn = header[i].split('/')[0]
                if cn == cname: i += 1
                else: break
            colls.append(CollTrimFunction(header[start:i], data[:, start:i]))

    return colls

def get_cfunctions(trimfile, fill, ignore=[]):
    CombCF = namedtuple("CombCF", "trim fill") # combined CollFunction

    fs = []
    tfuncs = cfunc_from_trimfile(trimfile)
    errors = 0
    ignored = 0
    for tf in tfuncs:
        if tf.cname in ignore:
            print("ignoring '{}'".format(tf.cname))
            ignored += 1
            continue

        cf = CombCF(tf, CollFillFunction(fill, tf.cname))
        try:
            cf.fill.fetch()
        except:
            print("error fetching fill")
            errors += 1
        else:
            fs.append(cf)
    print("---")
    print("{} collimator functions found".format(len(tfuncs)))
    print("\tignored: {}".format(ignored))
    print("\terrors:  {}".format(errors))
    print("\tin use:  {}".format(len(fs)))
    print("---")
    return fs
