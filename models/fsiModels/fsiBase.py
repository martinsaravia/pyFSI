from pyFSI.models.properties.dimensionlessNumbers import *
from abc import ABC, abstractmethod

class fsiBase(ABC):
    def __init__(self, execution, control, solid, flow, time):
        # ----- Public attributes ----- #
        self.base = "fsi"
        self.name = control['name']
        self.dof = None  # Number of degrees of freedom
        self.name = control['name']
        self.path = execution['paths']['fsiPath']  # Associated path
        self.output = []
        self.dimNumbers = {}  # Dimensionless numbers
        self.varMap = {
            "numbers": "dimNumbers"
        }

        # ----- Private attributes ----- #
        self._solid = solid
        self._flow = flow
        self._execution = execution
        self._control = control
        self._time = time
        if execution['debug'] == 'yes':
            self._debug = True
        else:
            self._debug = False

    def calcNumbers(self):
        self.dimNumbers["Cy"] = CauchyNumber(self)
        self.dimNumbers["Ms"] = massNumber(self)
        self.dimNumbers["Vr"] = reducedVelocityNumber(self)

    # Getters
    def control(self):
        return self._control

    def execution(self):
        return self._execution

    def solid(self):
        return self._solid

    def flow(self):
        return self._flow

    def time(self):
        return self._time

    def numbers(self):
        return self.dimNumbers

