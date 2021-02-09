from mesh.region.boundary1D import *

class boundary1DFlexible(boundary1D):
    def __repr__(self):
        return 'boundary1D: boundary1DFlexible'

    def __init__(self, mesh, name, y):
        super().__init__(mesh, name, coupled=True)

        # ----- Public attributes ----- #
        self.y0 = y  # Initial Position of the boundary
        self.y = y  # Initial Position of the boundary
        self.dy = np.zeros(mesh.size)  # Rigid boundaries have zero time derivative
        self.ddy = np.zeros(mesh.size)
        self.dyi = np.zeros(mesh.size)  # Indefinite integral of first time derivative
        self.ddyi = np.zeros(mesh.size)  # Indefinite integral of second time deriv
        self.vertices = None  # The 3D mesh

        self.mesh3D(mesh)

    def mesh3D(self, mesh):
        # Generate the 3D mesh of the interface
        self.vertices = np.zeros((mesh.size * 2, 3))  # The 3D mesh
        self.vertices[0:mesh.size, 0] = mesh.x
        self.vertices[0:mesh.size, 1] = self.y
        self.vertices[mesh.size:mesh.size*2, 0] = mesh.x
        self.vertices[mesh.size:mesh.size*2, 1] = self.y
        self.vertices[mesh.size:mesh.size*2, 2] = np.ones(mesh.size)  # This is the z=1 ghost nodes to math the solid 3d geometry
        # self.Forces = np.zeros_like(self.vertices)  # Force to transfer

    def setField(self, fieldName, data):
        if fieldName == "Displacements":
            self.y = self.y0 + data[0:self._mesh.size, 1]
            self.vertices[0:self._mesh.size, 1] = self.y
            self.vertices[self._mesh.size:self._mesh.size*2, 1] = self.y
        if fieldName == "Velocities":
            self.dy = data[0:self._mesh.size, 1]
            self.dyi = si.cumtrapz(self.dy, self._mesh.x, initial=0.0)
        if fieldName == "Accelerations":
            self.ddy = np.zeros_like(data[:, 0])
            self.ddyi = si.cumtrapz(self.ddy, self._mesh.x, initial=0.0)

    # ----- Abstract methods ----- #
    def update(self):
        # Re calculate the spatial derivatives and integrals
        self.calculate()

    def isFlexible(self):
        return True

    # ----- Getters ----- #
    def solid(self):
        return None

    def control(self):
        return self._control

    def eigen(self):
        return None
