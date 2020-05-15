# Change the system matrix.
# Corrected the position of G, Small acceleration and changed the sign of Tb

import numpy as np, scipy.integrate as si
from vectors.eigen import eigenSystem as es
import copy
from models.fsiModels.fsiModel import fsiModel

class lfb1D(fsiModel):
    def __repr__(self):
        return 'lfb1D '

    def __init__(self, execution, control, beam, flow):
        super().__init__(execution, control, beam, flow)

        # Time and parametric info
        self.ti = execution['time']['ti']
        self.pi = execution['parameters']['pi']

        # Size of the state-space eigen matrix
        self._size = beam.eigen()()['size'] * 2 + 2

        # Initialize the system matrices
        esize = beam.eigen()()['size']
        self.K = np.zeros((esize, esize))  # Mass
        self.C = np.zeros((esize, esize))  # Damping
        self.M = np.zeros((esize, esize))  # Stiffness
        self.Tb = np.zeros(esize)
        self.Tt = np.zeros(esize)
        self.Bb = np.zeros(esize)
        self.Bt = np.zeros(esize)
        self.Db = np.zeros(esize)
        self.Dt = np.zeros(esize)
        self.Eb = np.zeros(esize)
        self.Et = np.zeros(esize)
        self.Gt = None
        self.Gb = None
        self.norm = np.zeros((esize, esize))  # Eigenvector Norm
        self.S = np.zeros((self._size, self._size))  # System matrix
        self.ES = None

        # Assemble the fsi modal matrices K, C and M
        self.assemble(beam, flow)

        # Calculate the eigenvaluess and eigenvectors
        self.calculate()

    def assemble(self, beam, flow):
        # size of the eigensystem
        esize = beam.eigen()()['size']

        # Make the norm
        phiTphi = np.tensordot(beam.eigen().v(), beam.eigen().v(), axes=0)
        # Construct the norm matrix (assume orthogonality of modes, only diagonal
        # elements are calculated)
        for i in range(0, esize):
            for j in range(0, esize):
                self.norm[i, j] = si.simps(phiTphi[i, j], beam.mesh().x())

        # Fill the system Matrices
        for i in range(0, esize):
            gi = beam.eigen().v()[i]
            ki = beam.K()[i] + flow.K()[i]
            ci = beam.C()[i] + flow.C()[i]
            mi = beam.M()[i] + flow.M()[i]
            for j in range(0, esize):
                gj = beam.eigen().v()[j]
                self.K[i, j] = -si.simps(ki * gj, beam.mesh().x())
                self.C[i, j] = -si.simps(ci * gj, beam.mesh().x())
                self.M[i, j] = si.simps(mi * gj, beam.mesh().x())

        # System Region Vectors
        self.Gt = flow.Gq()['channelTop']
        self.Gb = flow.Gq()['channelBot']
        for i in range(0, esize):
            gi = beam.eigen().v()[i]
            self.Tt[i] = si.simps(flow.Tf()['channelTop'] * gi, beam.mesh().x())
            self.Tb[i] = si.simps(flow.Tf()['channelBot'] * gi, beam.mesh().x())  # Changed the sign
            self.Bt[i] = si.simps(flow.Bq()['channelTop'] * gi, beam.mesh().x())
            self.Bb[i] = si.simps(flow.Bq()['channelBot'] * gi, beam.mesh().x())
            self.Dt[i] = si.simps(flow.Dq()['channelTop'] * gi, beam.mesh().x())
            self.Db[i] = si.simps(flow.Dq()['channelBot'] * gi, beam.mesh().x())
            self.Et[i] = si.simps(flow.Eq()['channelTop'] * gi, beam.mesh().x())
            self.Eb[i] = si.simps(flow.Eq()['channelBot'] * gi, beam.mesh().x())

        # System  Matrix
        Mi = np.linalg.inv(self.M)

        self.S[0:esize, esize:2*esize] = np.identity(esize)
        self.S[esize:2*esize, 0:esize] = np.dot(Mi, self.K)
        self.S[esize:2*esize, esize:2*esize] = np.dot(Mi, self.C)
        self.S[esize:2*esize, 2*esize] = np.dot(Mi, self.Tb)
        self.S[esize:2*esize, 2*esize+1] = -np.dot(Mi, self.Tt)

        # Formulation considering acceleration terms
        if self._control['formulation'] == "Saravia":
            MiDotC = np.dot(Mi, self.C)
            MiDotK = np.dot(Mi, self.K)
            MiDotTb = np.dot(Mi, self.Tb)
            MiDotTt = np.dot(Mi, self.Tt)

            self.S[2 * esize, 0:esize] = -self.Eb - np.dot(self.Bb, MiDotK)
            self.S[2 * esize, esize:2 * esize] = -self.Db - np.dot(self.Bb, MiDotC)
            self.S[2 * esize, 2 * esize] = self.Gb - np.dot(self.Bb, MiDotTb)
            self.S[2 * esize, 2 * esize + 1] = np.dot(self.Bb, MiDotTt)

            self.S[2 * esize + 1, 0:esize] = self.Et + np.dot(self.Bt, MiDotK)
            self.S[2 * esize + 1, esize:2 * esize] = self.Dt + np.dot(self.Bt, MiDotC)
            self.S[2 * esize + 1, 2 * esize] = np.dot(self.Bt, MiDotTb)
            self.S[2 * esize + 1, 2 * esize + 1] = self.Gt - np.dot(self.Bt, MiDotTt)

        # Formulation not considering acceleration terms
        elif self._control['formulation'] == "SaraviaReduced":
            self.S[2 * esize, 0:esize] = -self.Eb
            self.S[2 * esize, esize:2 * esize] = -self.Db
            self.S[2 * esize, 2 * esize] = self.Gb

            self.S[2 * esize + 1, 0:esize] = self.Et
            self.S[2 * esize + 1, esize:2 * esize] = self.Dt
            self.S[2 * esize + 1, 2 * esize + 1] = self.Gt

        # Formulation of Tosi (I think some terms are wrong)
        elif self._control['formulation'] == "Tosi":
            S[2 * esize, 0:esize] = -(self.Eb + np.dot(self.Bb, np.dot(Mi, self.K)))
            S[2 * esize, esize:2 * esize] = -(self.Db + np.dot(self.Bb, np.dot(Mi, self.C)))
            S[2 * esize, 2 * esize] = -np.dot(self.Bb, np.dot(Mi, self.Tb))
            S[2 * esize, 2 * esize + 1] = self.Gb - np.dot(self.Bb, np.dot(Mi, self.Tt))

            S[2 * esize + 1, 0:esize] = self.Et + np.dot(self.Bt, np.dot(Mi, self.K))
            S[2 * esize + 1, esize:2 * esize] = self.Dt + np.dot(self.Bt, np.dot(Mi, self.C))
            S[2 * esize + 1, 2 * esize] = self.Gt + np.dot(self.Bt, np.dot(Mi, self.Tb))
            S[2 * esize + 1, 2 * esize + 1] = np.dot(self.Bt, np.dot(Mi, self.Tt))

    # Calculate eigenvalues and eigenvectors
    def calculate(self):
        # Calculate the eigenvalues and eigenvectors
        # ndof = solid.dict()['solution']['modes'] * 2 + 2
        evalues, evectors = np.linalg.eig(self.S)

        # Create an eigensystem object
        self.ES = es.eigenSystem(evalues, evectors, sort=True)




