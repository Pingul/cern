from enum import Enum


LogLevel = Enum("LogLevel", "info warning")

class ModuleLogger:

    def __init__(self, module_name):
        self.module_name = module_name

    def log(self, *args, log_level=LogLevel.info, end='\n', module_prestring=True):
        if log_level == LogLevel.warning:
            print("*WARNING*", end=" ")
        if module_prestring: 
            print(self.module_name, end=": ")
        print(*args, end=end)

if __name__ == "__main__":
    lg = ModuleLogger("Test")
    lg.log("testing message info", LogLevel.info)
    lg.log("testing message warning", LogLevel.warning)
    lg.log("testing default message")
    lg.log("testing", "many", "messages")

