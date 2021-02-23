import precice as prc

from pyFSI.models.fsiModels.fsiBase import fsiBase

class preCICE(fsiBase):
    def __repr__(self):
        return 'precice_FSI'

    def __init__(self, execution, control, solid, flow):
        super().__init__(execution, control, solid, flow)

        # Precice initialization
        print("Initializing the precice interface...")
        self.interface = prc.Interface("FLOW", "precice-config.xml", 0, 1)

        # Create the precice mesh objects (which are just pyFSI boundary objects)
        print("Creating the pyFSI-precice mesh objects...")
        self.meshes = []
        for mesh in self._control['meshes']:
            self.meshes.append(preciceMesh(mesh['name'],
                               mesh['read'],
                               mesh['write'],
                               self.interface,
                               self.flow().boundary()))

    def write(self):
        super.write()

    def finish(self):
        super.finish()

# A class to handle the mesh object in preCICE, actually it is a cloud of points
# on which the coupled variables are projected.
class preciceMesh:
    def __init__(self, name, read, write, interface, boundary):
        # Public Attributes
        self.name = name
        self.ID = interface.get_mesh_id(name)
        self.read = {}
        self.write = {}
        # Private attributes
        self._boundary = boundary[name]
        # Procedures
        # Build the dictionary of IDs for the fields
        for field in read:
            self.read[field] = interface.get_data_id(field, self.ID)
        for field in write:
            self.write[field] = interface.get_data_id(field, self.ID)

        # Set the 3D fluid mesh coordinates (the boundary coordinates made 3D)
        self.vertexIDs = interface.set_mesh_vertices(self.ID,
                                                     boundary[name].vertices)

    def readFields(self, interface):
        for fieldName, fieldID in self.read.items():
            data = interface.read_block_vector_data(fieldID,
                                                    self.vertexIDs)
            self._boundary.setField(fieldName,
                                    data)  # Store the field in the boundary object
        self._boundary.update()  # Update the boundary geometry and operators

    def writeFields(self, interface, flow):
        for fieldName, fieldID in self.write.items():
            if fieldName == "Forces":
                data = flow.Forces  # Pressure integral
                interface.write_block_vector_data(fieldID,
                                                  self.vertexIDs,
                                                  data)

    def boundary(self):
        return self._boundary