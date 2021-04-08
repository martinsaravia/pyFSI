# --------------------------------------------------------------------------- #
#    p    #     version: 0.2.4
#    y    #     date: 7/4/2021
#    F    #     author: Martin Saravia
#    S    #     description: Class for table reading
#    I    #     return: table object
# --------------------------------------------------------------------------- #
import numpy as np
from scipy import interpolate
from abc import ABC, abstractmethod


class fixedValueBase(ABC):
    def __init__(self, dictionary):
        self._dict = dictionary

    @abstractmethod
    def getValue(self, x=None):
        pass


# Fixed Value base class
class fixedValue(fixedValueBase):
    def __init__(self, dictionary):
        super().__init__(dictionary)

    def getValue(self, x=None):
        return self._dict['value']


# Linear interpolated boundary condition
class tabulatedFixedValue(fixedValueBase):
    def __init__(self, dictionary):
        super().__init__(dictionary)
        xAxis = self._dict["time"]
        yAxis = self._dict["value"]
        if "interpolation" in self._dict:
            kind = self._dict["interpolation"]
        else:
            kind = "linear"
        self._interpolator = interpolate.interp1d(xAxis, yAxis, fill_value='extrapolate', kind=kind)

    def getValue(self, x=None):
        return self._interpolator(x)

