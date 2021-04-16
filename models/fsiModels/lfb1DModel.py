# Change the system matrix.
# Corrected the position of G, Small acceleration and changed the sign of Tb
import copy
import numpy as np, scipy.integrate as si

from pyFSI.vectors.eigen import eigenSystem as es
from pyFSI.models.fsiModels.fsiBase import fsiBase
from pyFSI.models.properties.dimensionlessNumbers import dimensionlessNumber


class lfb1D(fsiBase):
    def __repr__(self):
        return 'lfb1D'

    def __init__(self, execution, control, solid, flow, time):
        super().__init__(execution, control, solid, flow, time)
        if self._debug:
            print("     WARNING: Region names are hardcoded.")

        # Size of the state-space eigen matrix
        esize = solid.eigen.size
        self._esize = esize
        self._size = self._esize * 2 + 2

        # Initialize the system matrices
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
        self.ES = es.eigenSystem([0], [0])  # Initialize the system to zero for correct file writing
        self.norm = np.zeros((esize, esize))  # Eigenvector Norm
        self.S = np.zeros((self._size, self._size))  # System matrix
        # Output variables
        self.varMap["eigenValues"] = "ES.evalues()"
        self.varMap["eigenVectors"] = "ES.evectors()"


    def update(self):
        # Update the flow, the solid is not updated
        self._flow.update()
        self.calcNumbers()
        self.assemble()

    def assemble(self):
        # Aliases
        solid = self._solid
        flow = self._flow
        esize = self._esize

        # Calculate the norm
        # phiTphi = np.tensordot(solid.eigen.vectors, solid.eigen.vectors, axes=0)
        # # Construct the norm matrix (assume orthogonality of modes, only diagonal
        # # elements are calculated)
        # for i in range(esize):
        #     for j in range(esize):
        #         self.norm[i, j] = si.simps(phiTphi[i, j], solid.mesh().x)

        # Fill the system Matrices
        for i in range(0, esize):
            gi = solid.eigen.vectors[i]
            ki = solid.k[i] + flow.k[i]
            ci = solid.c[i] + flow.c[i]
            mi = solid.m[i] + flow.m[i]

            # Galerkin discretization
            for j in range(0, esize):
                gj = solid.eigen.vectors[j]
                self.K[i, j] = -si.simps(ki * gj, solid.mesh().x)
                self.C[i, j] = -si.simps(ci * gj, solid.mesh().x)
                self.M[i, j] = si.simps(mi * gj, solid.mesh().x)

        # System Region Vectors
        self.Gt = flow.Gq['channelTop']
        self.Gb = flow.Gq['channelBot']

        # Galerkin discretization
        for i in range(0, esize):
            gi = solid.eigen.vectors[i]
            self.Tt[i] = si.simps(flow.Tf['channelTop'] * gi, solid.mesh().x)
            self.Tb[i] = si.simps(flow.Tf['channelBot'] * gi, solid.mesh().x)  # Changed the sign
            self.Bt[i] = si.simps(flow.Bq['channelTop'] * gi, solid.mesh().x)
            self.Bb[i] = si.simps(flow.Bq['channelBot'] * gi, solid.mesh().x)
            self.Dt[i] = si.simps(flow.Dq['channelTop'] * gi, solid.mesh().x)
            self.Db[i] = si.simps(flow.Dq['channelBot'] * gi, solid.mesh().x)
            self.Et[i] = si.simps(flow.Eq['channelTop'] * gi, solid.mesh().x)
            self.Eb[i] = si.simps(flow.Eq['channelBot'] * gi, solid.mesh().x)

        # Mass products
        Mi = np.linalg.inv(self.M)
        MiDotC = np.dot(Mi, self.C)
        MiDotK = np.dot(Mi, self.K)
        MiDotTb = np.dot(Mi, self.Tb)
        MiDotTt = np.dot(Mi, self.Tt)

        # System matrix
        self.S[0:esize, esize:2*esize] = np.identity(esize)
        self.S[esize:2*esize, 0:esize] = MiDotK
        self.S[esize:2*esize, esize:2*esize] = MiDotC
        self.S[esize:2*esize, 2*esize] = MiDotTb
        self.S[esize:2*esize, 2*esize+1] = -MiDotTt

        # Formulation considering acceleration terms
        if self._control['type'] == "Saravia":
            self.S[2 * esize, 0:esize] = -self.Eb - np.dot(self.Bb, MiDotK)
            self.S[2 * esize, esize:2 * esize] = -self.Db - np.dot(self.Bb, MiDotC)
            self.S[2 * esize, 2 * esize] = self.Gb - np.dot(self.Bb, MiDotTb)
            self.S[2 * esize, 2 * esize + 1] = np.dot(self.Bb, MiDotTt)

            self.S[2 * esize + 1, 0:esize] = self.Et + np.dot(self.Bt, MiDotK)
            self.S[2 * esize + 1, esize:2 * esize] = self.Dt + np.dot(self.Bt, MiDotC)
            self.S[2 * esize + 1, 2 * esize] = np.dot(self.Bt, MiDotTb)
            self.S[2 * esize + 1, 2 * esize + 1] = self.Gt - np.dot(self.Bt, MiDotTt)

        # Formulation not considering acceleration terms
        elif self._control['type'] == "SaraviaReduced":
            self.S[2 * esize, 0:esize] = -self.Eb
            self.S[2 * esize, esize:2 * esize] = -self.Db
            self.S[2 * esize, 2 * esize] = self.Gb

            self.S[2 * esize + 1, 0:esize] = self.Et
            self.S[2 * esize + 1, esize:2 * esize] = self.Dt
            self.S[2 * esize + 1, 2 * esize + 1] = self.Gt

        # Formulation of Tosi (I think some terms are wrong)
        elif self._control['type'] == "Tosi":
            self.S[2 * esize, 0:esize] = -(self.Eb + np.dot(self.Bb, np.dot(Mi, self.K)))
            self.S[2 * esize, esize:2 * esize] = -(self.Db + np.dot(self.Bb, np.dot(Mi, self.C)))
            self.S[2 * esize, 2 * esize] = -np.dot(self.Bb, np.dot(Mi, self.Tb))
            self.S[2 * esize, 2 * esize + 1] = self.Gb - np.dot(self.Bb, np.dot(Mi, self.Tt))

            self.S[2 * esize + 1, 0:esize] = self.Et + np.dot(self.Bt, np.dot(Mi, self.K))
            self.S[2 * esize + 1, esize:2 * esize] = self.Dt + np.dot(self.Bt, np.dot(Mi, self.C))
            self.S[2 * esize + 1, 2 * esize] = self.Gt + np.dot(self.Bt, np.dot(Mi, self.Tb))
            self.S[2 * esize + 1, 2 * esize + 1] = np.dot(self.Bt, np.dot(Mi, self.Tt))

        else:
            sys.exit("ERROR: No type in fsi formulation found...")

    def calcNumbers(self):
        super().calcNumbers()
        self.dimNumbers['Mr'] = massRatio(self)
        self.dimNumbers['Kr'] = stiffnessRatio(self)
        self.dimNumbers['Gr'] = gapRatio(self)
        self.dimNumbers['Vp'] = viscousParameter(self)


# Dimensional numbers of this model
class massRatio(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Mr"
        mf = fsi.flow().fluid()['rho'] * fsi.flow().lRef  # Fluid mass
        ms = fsi.solid().material()['rho'] * fsi.solid().tRef  # Solid mass
        self.value = ms / mf


class stiffnessRatio(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Kr"
        kf = fsi.flow().fluid()['rho'] * fsi.flow().Q0**2 * fsi.flow().lRef**3 / fsi.flow().dRef**2
        ks = fsi.solid().material()['E'] * fsi.solid().control()['I']   # Solid mass
        self.value = ks / kf


class gapRatio(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Gr"
        self.value = fsi.flow().eRef


class viscousParameter(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Vp"
        self.value = fsi.flow().eRef**2 * fsi.flow().dimNumbers["Re"].value
