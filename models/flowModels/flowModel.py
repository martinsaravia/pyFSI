from abc import ABCMeta, abstractmethod
from properties.materialProperties import fluids as db
from properties.dimensionlessNumbers import makeDimensionlessNumbers

# Base class for the fluid models

class flowModel(metaclass=ABCMeta):
    def __repr__(self):
        return 'fluidModel Abstract Class'

    def __init__(self, mesh, flowDict):
        # References to the mesh, boundary conditions and regions
        self._mesh = mesh
        self._bc = flowDict['bc']
        self._regions = flowDict['regions']  # Regions dict
        self._dict = flowDict

        # Get the fluid properties
        if "db" in flowDict["fluid"]:  # Specify the fluid from the database
            self._fluid = db[flowDict["fluid"]["db"]]
        else:
            self._fluid = flowDict["fluid"]  # Specify the properties directly

        # Initialize the reference parameters and the dimensionless numbers
        self.vRef = None  # Reference velocity
        self.lRef = None  # Reference length
        self.dRef = None  # Reference inlet size
        self.dimNumbers = None

    # Calculate the dimensionless numbers
    def calcNumbers(self):
        self.dimNumbers = makeDimensionlessNumbers(flow=self)

    # Getters
    def mesh(self):
        return self._mesh

    def dict(self):
        return self._dict

    def fluid(self):
        return self._fluid

    def bc(self):
        return self._bc

    def regions(self):
        return self._regions

