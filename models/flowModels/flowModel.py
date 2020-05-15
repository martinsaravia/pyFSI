from abc import ABCMeta, abstractmethod
from properties.materialProperties import fluids as db
from properties.dimensionlessNumbers import makeDimensionlessNumbers

# Base class for the fluid models

class flowModel(metaclass=ABCMeta):
    def __repr__(self):
        return 'fluidModel Abstract Class'

    def __init__(self, execution, control, mesh, boundary):
        # References to the mesh, boundary conditions and regions
        self._excution = execution
        self._control = control
        self._mesh = mesh
        self._boundary = boundary
        self._bc = control['bc']
        self._regions = control['regions']  # Regions dict

        # Get the fluid properties
        if "db" in control["fluid"]:  # Specify the fluid from the database
            self._fluid = db[control["fluid"]["db"]]
        else:
            self._fluid = control["fluid"]  # Specify the properties directly

        # Initialize the reference parameters and the dimensionless numbers
        self.vRef = None  # Reference velocity
        self.lRef = None  # Reference length
        self.dRef = None  # Reference inlet size
        self.dimNumbers = None

    # Calculate the dimensionless numbers
    def calcNumbers(self):
        self.dimNumbers = makeDimensionlessNumbers(flow=self)

    # Getters
    def execution(self):
        return self._execution

    def control(self):
        return self._control

    def mesh(self):
        return self._mesh

    def boundary(self):
        return self._boundary

    def fluid(self):
        return self._fluid

    def bc(self):
        return self._bc

    def regions(self):
        return self._regions

