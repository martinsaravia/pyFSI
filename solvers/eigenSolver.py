import numpy as np

from pyFSI.vectors.eigen import eigenSystem as es
from pyFSI.solvers.solverBase import solverBase

class eigen(solverBase):
    def __init__(self, fsi, odb):
        super().__init__(fsi, odb)

    def solve(self):
        time = self._time
        while time.value <= time.endTime:
            print("  Solving the MFSI case:", self._fsi.name, "for time ", time.value)
            self._fsi.update()
            evalues, evectors = np.linalg.eig(self._fsi.S)
            self._fsi.ES = es.eigenSystem(evalues, evectors, sort=True)
            self._odb.write()
            time.advance()

        self._odb.close()
