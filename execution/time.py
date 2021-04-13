import time as tt
import numpy as np

class time:
    def __init__(self, execution):
        # ----- Public attributes ----- #
        self.name = "time"
        self.path = execution['paths']['fsiPath']
        self.stepping = execution['time']['stepping']
        self.start = execution['time']['startTime']
        self.end = execution['time']['endTime']
        self.delta = execution['time']['deltaT']
        self.value = self.start  # Current time
        self.startDate = tt.perf_counter()
        self.span = np.array([self.start, self.start + self.delta])
        # Output variable mapping
        self.varMap = {
            "time":       "value"
        }

    def elapsed(self):
        return tt.perf_counter() - self.startDate

    def advance(self):
        self.value += self.delta
        self.span += self.delta
