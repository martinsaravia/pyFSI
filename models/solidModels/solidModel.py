from abc import ABCMeta
from properties.materialProperties import solids as dataBase
from properties.dimensionlessNumbers import makeDimensionlessNumbers

class solidModel(metaclass=ABCMeta):
    def __repr__(self):
        return 'solidModel Abstract Class'

    def __init__(self, execution, control, mesh, name=None):
        # References to dictionaries and mesh
        self._execution = execution
        self._control = control
        self._mesh = mesh
        self.name = name

        # Get the material properties
        if 'db' in control['material']:  # Specify the material from the database
            self._material = dataBase[control['material']['db']]
        else:
            self._material = control['material']  # Specify the properties directly

        # Initialize variables
        self.lRef = None  # Reference length
        self.uRef = None  # Reference displacement
        self.vRef = None  # Reference velocity
        self.dimNumbers = None  # Dimensionless numbers

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
