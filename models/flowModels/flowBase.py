"""@package docstring
Base class for the fluid models

Only one boundary can be flexible
"""

from abc import ABCMeta, abstractmethod
from pyFSI.models.properties.materialProperties import fluids as db
from pyFSI.models.properties.dimensionlessNumbers import makeDimensionlessNumbers
from pyFSI.mesh.region.fsiRegion1D import fsiRegion1D


class flowModel(metaclass=ABCMeta):
    """ Base class for the flow models"""
    def __repr__(self):
        return 'fluidModel Abstract Class'

    def __init__(self, execution, control, mesh, boundary):
        """Constructor"""

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
        if execution['debug'] == 'yes':
            self._debug = True
        else:
            self._debug = False

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
        """ Calculates dimensionless numbers"""
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

    def closeOutput(self):
        for i in self.output:
            i.close()  # Close all files

    def finish(self):
        self.closeOutput()

