from abc import ABCMeta
from pyFSI.models.properties.materialProperties import solids as dataBase
from pyFSI.models.properties.dimensionlessNumbers import makeDimensionlessNumbers

class solidModel(metaclass=ABCMeta):
    def __repr__(self):
        return 'solidModel Abstract Class'

    def __init__(self, execution, control, mesh, name=None):

        # ----- Public attribues ----- #
        self.name = name
        self.dof = None
        self.output = []  # List of file objects
        self.lRef = None  # Reference length
        self.uRef = None  # Reference displacement
        self.vRef = None  # Reference velocity
        self.dimNumbers = None  # Dimensionless numbers

        # ----- Private attributes ----- #
        self._execution = execution
        self._control = control
        self._mesh = mesh

        # ----- Procedures ----- #
        # Get the material properties
        if 'db' in control['material']:  # Specify the material from the database
            self._material = dataBase[control['material']['db']]
        else:
            self._material = control['material']  # Specify the properties directly

    # Calculate the dimensionless numbers
    def calcNumbers(self):
        self.dimNumbers = makeDimensionlessNumbers(solid=self)

    # Getters
    def execution(self):
        return self._execution

    def control(self):
        return self._control

    def mesh(self):
        return self._mesh

    def material(self):
        return self._material

    def finish(self):
        for i in self.output:
            i.close()  # Close all files
