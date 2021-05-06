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
        self._execution = fsi.execution()

        # Public Attributes
        self.control = fsi.execution()['solver']

        # Output
        bufferSize = 1
        self.output = []
        self.output.append(open(self._execution['paths']['fsiPath'] / 'time.out', 'a+', buffering=bufferSize))

    # Abstract methods
    @abstractmethod
    def solve(self):
        pass

    # Getters
    def execution(self):
        return self._execution

    def odb(self):
        return self._odb

