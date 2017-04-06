from enum import Enum


LogLevel = Enum("LogLevel", "info warning")

class ModuleLogger:

    def __init__(self, module_name):
        self.module_name = module_name

    def log(self, message, log_level=LogLevel.info):
        if log_level == LogLevel.warning:
            print("*WARNING*", end=" ")
        print("{}: {}".format(self.module_name, message))

if __name__ == "__main__":
    lg = ModuleLogger("Test")
    lg.log("testing message info", LogLevel.info)
    lg.log("testing message warning", LogLevel.warning)
    lg.log("testing default message")

