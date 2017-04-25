# encoding: utf-8
from enum import Enum

LogLevel = Enum("LogLevel", "info success notify warning")

def color_wrap(string, color):
    return '\x1b[{}m'.format(color) + string + '\x1b[0m'

class ModuleLogger:
    colors = {
            LogLevel.notify : '1;36;40',
            LogLevel.warning : '1;31;40',
            LogLevel.success : '1;32;40'
            }

    def __init__(self, module_name):
        self.module_name = module_name

    def log(self, *args, log_level=LogLevel.info, end='\n', module_prestring=True):
        message = " ".join(map(str, [*args]))
        if module_prestring:
            message = "{}: {}".format(self.module_name, message)
            if log_level == LogLevel.warning:
                message = "*WARNING* {}".format(message)
        
        if log_level in self.colors:
            message = color_wrap(message, self.colors[log_level])

        print(message, end=end)

if __name__ == "__main__":
    # Some testcases
    lg = ModuleLogger("Test")
    lg.log("message info", log_level=LogLevel.info)
    lg.log("message success", log_level=LogLevel.success)
    lg.log("message notify", log_level=LogLevel.notify)
    lg.log("message warning", log_level=LogLevel.warning)
    lg.log("testing default message")
    lg.log("testing", "many", "messages")

    lg.log("message", end=' ')
    lg.log("that", end=' ', module_prestring=False)
    lg.log("continues", module_prestring=False)
