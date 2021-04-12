import time as tt


class time:
    def __init__(self, execution):
        # ----- Public attributes ----- #
        self.name = "time"
        self.path = execution['paths']['fsiPath']
        self.stepping = execution['time']['stepping']
        self.startTime = execution['time']['startTime']
        self.endTime = execution['time']['endTime']
        self.deltaT = execution['time']['deltaT']
        self.value = self.startTime
        self.startDate = tt.perf_counter()
        # Output variable mapping
        self.varMap = {
            "value":       "value"
        }

    def elapsed(self):
        return tt.perf_counter() - self.startDate

    def advance(self):
        self.value += self.deltaT
