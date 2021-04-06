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
    def __init__(self, fsi):
        # Public Attributes
        self.fsi = fsi
        self.time = None
        self.control = fsi.execution()['solver']
        self._execution = fsi.execution()

    # Execute final tasks for the fsi objects
    def finish(self):
        self.fsi.finish()

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

