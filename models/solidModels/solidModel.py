from abc import ABCMeta
from properties.materialProperties import solids as dataBase
from properties.dimensionlessNumbers import makeDimensionlessNumbers

class solidModel(metaclass=ABCMeta):
    def __repr__(self):
        return 'solidModel Abstract Class'

    def __init__(self, mesh, controlDict, solidDict, name=None):
        # References to the mesh, control, and dictionary
        self._mesh = mesh
        self._control = controlDict
        self._dict = solidDict
        self.name = name

        # Get the material properties
        if 'db' in solidDict['material']:  # Specify the material from the database
            self._material = dataBase[solidDict['material']['db']]
        else:
            self._material = solidDict['material']  # Specify the properties directly

        # Initialize variables
        self.lRef = None  # Reference length
        self.uRef = None  # Reference displacement
        self.vRef = None  # Reference velocity
        self.dimNumbers = None  # Dimensionless numbers

    # Calculate the dimensionless numbers
    def calcNumbers(self):
        self.dimNumbers = makeDimensionlessNumbers(solid=self)

    # Getters
    def mesh(self):
        return self._mesh

    def dict(self):
        return self._dict

    def material(self):
        return self._material
