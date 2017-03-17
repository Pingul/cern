import json

# DIR = "calc/"
# # DIR = "simulations/store/5k_H8000/"

# PARTICLE_FILE = DIR + "particles.dat"
# LINE_FILE = DIR + "lines.dat"
# COLL_FILE = DIR + "coll.dat"
# STARTDIST_FILE = DIR + "startdist.dat"
# ENDDIST_FILE = DIR + "enddist.dat"
# RAMP_FILE = "resources/LHC_ramp.dat"

# ## 1 bucket
# PLOT_FRAME = {
    # 'x' : [-0.1*pi, 2.1*pi],
    # 'y' : [-0.4e9, 0.4e9]
# }

# ## 3 buckets
# # PLOT_FRAME = {
    # # 'x' : [-2*pi, 4*pi],
    # # 'y' : [-2e9, 2e9]
# # }

# ## Constants
# H_SEPARATRIX = 1909859.317103239 # Taken from 2d-synchrotron.hh
# F_RF = 400788731.3727857
# BLM_INT = int(1.3*11245.0)

SETTINGS_FILE = "settings.json"
class Settings(dict):
    """ dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __presave(self):
        del self.PARTICLE_PATH
        del self.LINE_PATH 
        del self.COLL_PATH 
        del self.STARTDIST_PATH
        del self.ENDDIST_PATH
        del self.RAMP_PATH 

    def __postload(self):
        # Special cases
        self.PARTICLE_PATH = "{}/{}".format(self.DATA_DIR, self.PARTICLE_FILE)
        self.LINE_PATH = "{}/{}".format(self.DATA_DIR, self.LINE_FILE)
        self.COLL_PATH = "{}/{}".format(self.DATA_DIR, self.COLL_FILE)
        self.STARTDIST_PATH = "{}/{}".format(self.DATA_DIR, self.STARTDIST_FILE)
        self.ENDDIST_PATH = "{}/{}".format(self.DATA_DIR, self.ENDDIST_FILE)
        self.RAMP_PATH = "{}/{}".format(self.RESOURCE_DIR, self.RAMP_FILE)

    def save(self):
        self.presave()
        print("reading settings from '{}'".format(SETTINGS_FILE))
        with open(SETTINGS_FILE, 'w') as f:
            json.dump(self, f)
    
    @classmethod
    def load(clss):
        print("loading settings from '{}'".format(SETTINGS_FILE))
        with open(SETTINGS_FILE, 'r') as f:
            s = Settings(json.load(f))
            s.__postload()
            return s

settings = Settings.load()
