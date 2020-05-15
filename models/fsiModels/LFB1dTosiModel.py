import numpy as np, scipy.integrate as si
from vectors.eigen import eigenSystem as es
from models.fsiModels.fsiModel import fsiModel

class LFB1dTosi(fsiModel):
    def __repr__(self):
        return 'LFB1dTosi'

    def __init__(self, control, beam, fluid):

        # Time and parametric info
        self.ti = control['time']['ti']
        self.pi = control['parameters']['pi']

        self.beam = beam
        self.fluid = fluid
        # Size of the state-space eigen matrix
        self._size = beam.eigen()()['size'] * 2 + 2

        # Store the control dict
        self.control = control

        # Assemble the fsi modal matrices K, C and M
        self.assemble(beam, fluid)

        # Calculate the eigenvaluess and eigenvectors
        self.calculate()


    def assemble(self, beam, fluid):
        # size of the eigensystem
        esize = beam.eigen()()['size']

        # Make the norm
        self.norm = np.zeros((esize, esize)) # Eigenvector Norm
        phiTphi = np.tensordot(beam.eigen().v(), beam.eigen().v(), axes=0)

        # Construct the norm matrix (assume orthogonality of modes, only diagonal
        # elements are calculated)
        for i in range(0, esize):
            for j in range(0, esize):
                self.norm[i, j] = si.simps(phiTphi[i, j], beam.mesh().x())

        # System Matrices
        self.K = np.zeros((esize, esize))  # Mass
        self.C = np.zeros((esize, esize))  # Mass
        self.M = np.zeros((esize, esize))  # Mass
        for i in range(0, esize):
            gi = beam.eigen().v()[i]
            ki = beam.K()[i] + fluid.K()[i]
            ci = beam.C()[i] + fluid.C()[i]
            mi = beam.M()[i] + fluid.M()[i]
            for j in range(0, esize):
                gj = beam.eigen().v()[j]
                self.K[i, j] = -si.simps(ki * gj, beam.mesh().x())
                self.C[i, j] = -si.simps(ci * gj, beam.mesh().x())
                self.M[i, j] = si.simps(mi * gj, beam.mesh().x())

        # System Region Vectors
        self.Tb = np.zeros(esize)
        self.Tt = np.zeros(esize)
        self.Bb = np.zeros(esize)
        self.Bt = np.zeros(esize)
        self.Db = np.zeros(esize)
        self.Dt = np.zeros(esize)
        self.Eb = np.zeros(esize)
        self.Et = np.zeros(esize)
        self.Gt = fluid.Gq()['channelTop']
        self.Gb = fluid.Gq()['channelBot']
        for i in range(0, esize):
            gi = beam.eigen().v()[i]
            self.Tt[i] = -si.simps(fluid.Tf()['channelTop'] * gi, beam.mesh().x())
            self.Tb[i] = si.simps(fluid.Tf()['channelBot'] * gi, beam.mesh().x())

            self.Bt[i] = si.simps(fluid.Bq()['channelTop'] * gi, beam.mesh().x())
            self.Bb[i] = si.simps(fluid.Bq()['channelBot'] * gi, beam.mesh().x())

            self.Dt[i] = si.simps(fluid.Dq()['channelTop'] * gi, beam.mesh().x())
            self.Db[i] = si.simps(fluid.Dq()['channelBot'] * gi, beam.mesh().x())
            self.Et[i] = si.simps(fluid.Eq()['channelTop'] * gi, beam.mesh().x())
            self.Eb[i] = si.simps(fluid.Eq()['channelBot'] * gi, beam.mesh().x())


        # System  Matrix
        Mi = np.linalg.inv(self.M)
        gsize = esize * 2 + 2
        S = np.zeros((gsize, gsize))

        S[0:esize, esize:2*esize] = np.identity(esize)
        S[esize:2*esize, 0:esize] = np.dot(Mi, self.K)
        S[esize:2*esize, esize:2*esize] = np.dot(Mi, self.C)
        S[esize:2*esize, 2*esize] = np.dot(Mi, self.Tb)
        S[esize:2*esize, 2*esize+1] = np.dot(Mi, self.Tt)

        S[2*esize, 0:esize] = -(self.Eb + np.dot(self.Bb, np.dot(Mi, self.K)))
        S[2*esize, esize:2*esize] = -(self.Db + np.dot(self.Bb, np.dot(Mi, self.C)))
        S[2*esize, 2*esize] = -np.dot(self.Bb, np.dot(Mi, self.Tb))
        S[2*esize, 2*esize+1] = self.Gb - np.dot(self.Bb, np.dot(Mi, self.Tt))

        S[2*esize+1, 0:esize] = self.Et + np.dot(self.Bt, np.dot(Mi, self.K))
        S[2*esize+1, esize:2*esize] = self.Dt + np.dot(self.Bt, np.dot(Mi, self.C))
        S[2*esize+1, 2*esize] = self.Gt + np.dot(self.Bt, np.dot(Mi, self.Tb))
        S[2*esize+1, 2*esize+1] = np.dot(self.Bt, np.dot(Mi, self.Tt))

        self.S = S

        # Version without loop and not normalized
        # mm = beam.mass() + fluid.mass() # Unintegrated Mass
        # self.M = si.simps(np.tensordot(mm.T, es.v(), axes=0))#2.97 LFB1dTosi

    # Calculate eigenvalues and eigenvectors
    def calculate(self):
        # Calculate the eigenvalues and eigenvectors
        # ndof = solid.dict()['solution']['modes'] * 2 + 2
        evalues, evectors = np.linalg.eig(self.S)

        # Create an eigensystem object
        self._ES  = es.eigenSystem(evalues, evectors, sort=True)

    # Get the eigensystem
    def eigen(self):
        return self._ES
