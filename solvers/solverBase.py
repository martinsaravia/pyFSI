# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 07/12/2020
#    F    #     author: Martin Saravia
#    S    #     description: Base class for all solvers
#    I    #
# --------------------------------------------------------------------------- #

import importlib
from abc import ABCMeta, abstractmethod

from pyFSI.mesh.fsiMesh1D import fsiMesh1D


class solverBase:
    def __init__(self, fsi, odb):
        # Private Attributes
        self._fsi = fsi
        self._odb = odb
        self._time = fsi.time()
        self._control = fsi.execution()['solver']
        self._execution = fsi.execution()

        # Output
        bufferSize = 1
        self.output = []
        self.output.append(open(self._execution['paths']['fsiPath'] / 'time.out', 'a+', buffering=bufferSize))

    # Abstract methods
    @abstractmethod
    def solve(self):
        pass

    @abstractmethod
    def advance(self, tspan):
        pass

    # Write the output of every object
    @abstractmethod
    def write(self):
        pass

    # Getters
    def execution(self):
        return self._execution

    def odb(self):
        return self._odb

