import json

class SettingsBase(dict):
    """ dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    
    def presave(self):
        pass

    def postload(self):
        pass

    @classmethod
    def settings_file(clss):
        raise Exception("Needs to be implemented")


    @classmethod
    def save(clss, settings):
        s = clss(settings)
        s.presave()
        print("writing settings to '{}'".format(clss.settings_file()))
        with open(clss.settings_file(), 'w') as f:
            json.dump(s, f, indent=4, separators=(',', ': '), sort_keys=True)
    
    @classmethod
    def load(clss):
        print("reading settings from '{}'".format(clss.settings_file()))
        with open(clss.settings_file(), 'r') as f:
            s = clss(json.load(f))
            s.postload()
            return s
