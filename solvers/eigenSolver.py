import numpy as np

from pyFSI.vectors.eigen import eigenSystem as es
from pyFSI.solvers.solverBase import solverBase

class eigen(solverBase):
    def __init__(self, fsi):
        super().__init__(fsi)

    def solve(self):
        print("  Solving ", self.fsi.name)
        evalues, evectors = np.linalg.eig(self.fsi.S)
        self.solution = es.eigenSystem(evalues, evectors, sort=True)
        self.write(self.solution)

    def write(self, solution):
        self.fsi.write(solution)
        self.fsi.finish()


