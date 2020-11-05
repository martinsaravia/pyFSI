# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: BE Beam class
#    I    #     return: solidModel object
# --------------------------------------------------------------------------- #
# Notes:
#   Nonlinear FSI beam-lf model
#
# --------------------------------------------------------------------------- #
# Change the system matrix.
# Corrected the position of G, Small acceleration and changed the sign of Tb

import numpy as np, scipy.integrate as si
from vectors.eigen import eigenSystem as es
import copy
from models.fsiModels.fsiModel import fsiModel
from properties.dimensionlessNumbers import dimensionlessNumber
import pandas as pd

class nlLFB1D(fsiModel):
    def __repr__(self):
        return 'nlLFB1D '

    def __init__(self, execution, control, beam, flow):
        super().__init__(execution, control, beam, flow)

        # ----- Public attributes ----- #
        self.dof = 2 * beam.dof + flow.dof  # DOFs of the FSI system
        self.state = None  # Current State
        self.rhs = np.zeros(self.dof)  # Current variables velocities
        self.sol = []
        self.out = {'t': [],
                    'state': []}

        # ----- Private attributes ----- #
        # Indexing coefficients
        self._i0 = 0
        self._i1 = beam.dof
        self._i2 = 2 * beam.dof
        self._i3 = self.dof

        # ----- Procedures ----- #
        self.setInitialConditions(beam, flow)

    def setInitialConditions(self, beam, flow):
        beam.setInitialConditions()
        flow.setInitialConditions()
        # Inialize the state variables
        state0 = np.zeros(self.dof)  # Current initial variables values
        state0[self._i0:self._i1] = beam.a
        state0[self._i1:self._i2] = beam.da
        state0[self._i2:self._i3] = flow.Q0
        self.state = state0
        self.sState = state0[self._i0:self._i2] # Solid state
        self.fState = state0[self._i2:self._i3] # Fluid State

    def update(self, who, time, state):
        # self.sol.append(state)
        if who == 'solid':
            self.sState = state
            # Update the solid (UPDATED FIRST BECAUSE IT IS A FLUID BOUNDARY)
            a = state[self._i0:self._i1]
            da = state[self._i1:self._i2]
            # dda = self.rhsSolid(None, state)[self._i1:self._i2] # Get the acceleration BAD????
            self._beam.update(a, da, self.rhs[self._i0:self._i2])
            # Write the results to file
            self._beam.write()

        elif who == 'fluid':
            # Update the flow
            # Update the rhs (this is, dQ) with the new Q and get the pressure
            # right since it is a function of dQ
            self.fState = state
            # self._flow.rhs(time, state)
            self._flow.update(time, state)
            self._flow.write()

        elif who == 'bc':
            self._flow.pUpdate(time)

        elif who == 'regions':
            # Update the regions with the new beam deformation
            for i, region in enumerate(self._flow.regions):
                region.update()

    def rhsSolid(self, time, solidState):
        beam = self._beam
        flow = self._flow

        self.rhs[self._i0:self._i2] = np.dot(beam.S, solidState) + beam.F + beam.addedStateModalForce(flow.deltaPx[0])

        # print('beamF', beam.F)
        # print('flowF', beam.addedStateModalForce(flow.deltaPx[0]))

        return self.rhs[self._i0:self._i2]


    # def residual(self):
    #     res = self._beam.F() + self._flow.F() - self.beam.

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


