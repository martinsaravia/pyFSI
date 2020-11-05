# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: Flow model interface
#    I    #
# --------------------------------------------------------------------------- #
# Notes:
#   This class generates a region from two 1D boundary objects
#   Only one boundary can be flexible
#
# --------------------------------------------------------------------------- #
from abc import ABCMeta, abstractmethod
from properties.materialProperties import fluids as db
from properties.dimensionlessNumbers import makeDimensionlessNumbers
from mesh.region.fsiRegion1D import fsiRegion1D
# Base class for the fluid models

class flowModel(metaclass=ABCMeta):
    def __repr__(self):
        return 'fluidModel Abstract Class'

    def __init__(self, execution, control, mesh, boundary):
        # ----- Public attribues ----- #
        self.dof = None  # Number of DOF of the model
        self.regions = []  # Regions comprising the domain
        self.vRef = None  # Reference velocity
        self.lRef = None  # Reference length
        self.dRef = None  # Reference inlet size
        self.dimNumbers = None
        self.output = []

        # ----- Private attributes ----- #
        self._execution = execution
        self._control = control
        self._mesh = mesh
        self._boundary = boundary
        self._fluid = None

        # ----- Procedures ----- #
        # Get the fluid properties
        if "db" in control["fluid"]:  # Specify the fluid from the database
            self._fluid = db[control["fluid"]["db"]]
        else:
            self._fluid = control["fluid"]  # Specify the properties directly

        # Create the region objects
        [self.regions.append(fsiRegion1D(i, mesh, boundary)) for i in control['regions']]

    # Calculate the dimensionless numbers
    def calcNumbers(self):
        self.dimNumbers = makeDimensionlessNumbers(flow=self)

    # Pure virtual methods

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

    def finish(self):
        for i in self.output:
            i.close()  # Close all files

