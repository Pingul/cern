import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../common/"))
from settings_base import SettingsBase

class Settings(SettingsBase):

    def presave(self):
        del self.PARTICLE_PATH
        del self.LINE_PATH 
        del self.COLL_PATH 
        del self.STARTDIST_PATH
        del self.ENDDIST_PATH
        del self.RAMP_PATH 
        del self.META_PATH

    def postload(self):
        # Special cases
        self.META_PATH = "{}/{}".format(self.DATA_DIR, self.META_FILE)
        self.PARTICLE_PATH = "{}/{}".format(self.DATA_DIR, self.PARTICLE_FILE)
        self.LINE_PATH = "{}/{}".format(self.DATA_DIR, self.LINE_FILE)
        self.COLL_PATH = "{}/{}".format(self.DATA_DIR, self.COLL_FILE)
        self.STARTDIST_PATH = "{}/{}".format(self.DATA_DIR, self.STARTDIST_FILE)
        self.ENDDIST_PATH = "{}/{}".format(self.DATA_DIR, self.ENDDIST_FILE)
        self.RAMP_PATH = "{}/{}".format(self.RESOURCE_DIR, self.RAMP_FILE)

    @classmethod
    def settings_file(clss):
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "settings.json")

settings = Settings.load()

