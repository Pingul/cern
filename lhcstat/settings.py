import sys, os
sys.path.append("../common/")
from settings_base import SettingsBase

class OMLSettings(SettingsBase):
    @classmethod
    def settings_file(clss):
        return os.path.join(os.path.dirname(__file__), 'settings.json')

settings = OMLSettings.load()
