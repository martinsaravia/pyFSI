# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: BE solid class
#    I    #     return: solidModel object
# --------------------------------------------------------------------------- #
# Notes:
#   Nonlinear FSI solid-lf model
#
# --------------------------------------------------------------------------- #
# Change the system matrix.
# Corrected the position of G, Small acceleration and changed the sign of Tb
import copy
import numpy as np
import pandas as pd
import scipy.integrate as si

from pyFSI.models.fsiModels.fsiBase import fsiBase
from pyFSI.models.properties.dimensionlessNumbers import dimensionlessNumber
from pyFSI.vectors.eigen import eigenSystem as es


class nlLFB1D(fsiBase):
    def __repr__(self):
        return 'nlLFB1D '

    def __init__(self, execution, control, solid, flow):
        super().__init__(execution, control, solid, flow)

        # ----- Public attributes ----- #
        self.dof = 2 * solid.dof + flow.dof  # DOFs of the FSI system
        self.state = None  # Current State
        self.fState = None
        self.sState = None
        self.rhs = np.zeros(self.dof)  # Current variables velocities

        # ----- Private attributes ----- #
        # Indexing coefficients
        self._i0 = 0
        self._i1 = solid.dof
        self._i2 = 2 * solid.dof
        self._i3 = self.dof

        # ----- Procedures ----- #
        self.setInitialConditions(solid, flow)

    def setInitialConditions(self, solid, flow):
        solid.setInitialConditions()
        flow.setInitialConditions()
        # Inialize the state variables
        state0 = np.zeros(self.dof)  # Current initial variables values
        state0[self._i0:self._i1] = solid.a
        state0[self._i1:self._i2] = solid.da
        state0[self._i2:self._i3] = flow.Q0
        self.state = state0
        self.sState = state0[self._i0:self._i2]  # Solid state
        self.fState = state0[self._i2:self._i3]  # Fluid State

    def update(self, who, time, state):
        if who == 'solid':
            self.sState = state
            # Update the solid (UPDATED FIRST BECAUSE IT IS A FLUID BOUNDARY)
            a = state[self._i0:self._i1]
            da = state[self._i1:self._i2]
            # dda = self.rhsSolid(None, state)[self._i1:self._i2] # Get the acceleration WRONG????
            self._solid.update(a, da, self.rhs[self._i0:self._i2])
            # Write the results to file
            # self._solid.write()

        elif who == 'flow':
            # Update the flow
            # Update the rhs (this is, dQ) with the new Q and get the pressure
            # right since it is a function of dQ
            self.fState = state
            # self._flow.rhs(time, state)
            self._flow.update(time, state)  # Update the flow variables

            self._flow.updateForces(time)   # Update the

        elif who == 'regions':
            # Update the regions with the new solid deformation
            for i, region in enumerate(self._flow.regions):
                region.update()

    def rhsSolid(self, time, solidState):
        solid = self._solid
        flow = self._flow
        self.rhs[self._i0:self._i2] = (np.dot(solid.S, solidState) +
                                       solid.F +
                                       solid.addedStateModalForce(flow.deltaPx[0]))
        return self.rhs[self._i0:self._i2]

    # def residual(self):
    #     res = self._solid.F() + self._flow.F() - self.solid.

    def calcNumbers(self):
        super().calcNumbers()
        self.dimNumbers['Mr'] = massRatio(self)
        self.dimNumbers['Kr'] = stiffnessRatio(self)
        self.dimNumbers['Gr'] = gapRatio(self)
        self.dimNumbers['Vp'] = viscousParameter(self)

    def write(self):
        self._solid.write()
        self._flow.write()

    def finish(self):
        super().finish()


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
        kf = fsi.flow().fluid()['rho'] * fsi.flow().qx0**2 * fsi.flow().lRef**3 / fsi.flow().dRef**2
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
