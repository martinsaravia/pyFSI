# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: NL Leakage Flow Class
#    I    #     return: flow object
# --------------------------------------------------------------------------- #
# Notes:
#   Nonlinear leakage flow model
#
#
# Warnings:
#   Local Re is the same for all regions. Not as above Eq. 2.30
#   eta calculation must be checked
#   Check the consideration of the beam thickness (hardcoded to 1)
#   Only two regions allowed
#   Regions has just one flexible boundary
#   Inconsistency: self._control['bc']['inlet']['zeta'] is the same for regions
# --------------------------------------------------------------------------- #

from abc import ABC
import numpy as np
import scipy.integrate as si

from pyFSI.mesh.region.fsiRegion1D import fsiRegion1D
from pyFSI.models.flowModels.flowBase import flowModel
from pyFSI.models.properties.boundaryLayer import boundaryLayer
from pyFSI.vectors.eigen.eigenVector import eigenVector
from pyFSI.models.properties.dimensionlessNumbers import ReynoldsNumber

class nlLeakageFlow2D(flowModel, ABC):
    def __repr__(self):
        return 'leakageFlow2DModel '

    def __init__(self, execution, control, mesh, boundary,  name='NN'):
        super().__init__(execution, control, mesh, boundary)

        # ----- Public attributes ----- #
        self.dof = 2  # Number of Q equations (1 per region)
        self.Q0 = np.zeros(self.dof)
        self.v0 = np.zeros(self.dof)
        self.dQ0 = np.zeros(self.dof)
        self.Q = np.zeros(self.dof, dtype=object)
        self.dQ = np.zeros(self.dof, dtype=object)
        self.px = np.zeros(self.dof, dtype=object)
        self.deltaPx = np.zeros(1, dtype=object)
        self.Forces = np.zeros((mesh.size*2, 3))  # The 3D Forces
        self.Dp = None  # pOut - pIn
        self.converged = [False] * self.dof
        # self.eRef = np.zeros(self.dof)
        # self.bLayer = np.zeros(self.dof)  # Boundary layer

        # ----- Private attributes ----- #
        self._ti = None  # Current time
        self._th = 1  # Thickness of the model
        self._pIn = None  # Inlet pressure
        self._pOut = None  # Outlet pressure
        self._pTol = 1E-4  # Pressure convergence tolerance
        self._zetaIn = control['bc']['inlet']['zeta']  # Inlet loss factor
        self._zetaOut = control['bc']['outlet']['zeta']  # Outlet loss factor
        # self._f0 = np.zeros(self.dof)  # Viscous friction factor
        # self._xix = np.zeros(self.dof)  # Nonlinear profile factor
        # self._eta = np.zeros(self.dof)  # Derivative of f(Qx) at qx0
        self.dRef = np.zeros(self.dof)  # Reference width
        self.lRef = mesh.L  # Fixed Reference length
        self.vRef = np.zeros(self.dof)  # Reference velocity
        self.eRef = np.zeros(self.dof)  # Reference something
        self.f0 = np.zeros(self.dof)
        self.xix = np.zeros(self.dof)
        self.eta = np.zeros(self.dof)
        self.Rd = [None] * self.dof  # Emtpy list of Reynolds objects
        # Size of the associated eigensystem
        # self._gDof = self.regions[0].eigen().size

        # Output
        bufferSize = 1
        self.output = []
        self.output.append(open(self._execution['paths']['fluidPath'] / 'Q0.out',  'a+', buffering=bufferSize))
        self.output.append(open(self._execution['paths']['fluidPath'] / 'dQ0.out', 'a+', buffering=bufferSize))
        self.output.append(open(self._execution['paths']['fluidPath'] / 'p0.out',  'a+', buffering=bufferSize))
        self.output.append(open(self._execution['paths']['fluidPath'] / 'pIn.out', 'a+', buffering=bufferSize))
        self.output.append(open(self._execution['paths']['fluidPath'] / 'v0.out',  'a+', buffering=bufferSize))
        self.output.append(open(self._execution['paths']['fluidPath'] / 'time.out', 'a+',buffering=bufferSize))
        self.output.append(open(self._execution['paths']['fluidPath'] / 'force.out', 'a+', buffering=bufferSize))


    # Flow rate equation initial condition
    def setInitialConditions(self):
        print("---> Starting Q0 iteration with initial Q0 = ", self.Q0)
        # Initialize the pressure difference
        self._pOut = self._control['bc']['outlet']['p']
        self._pIn = self._control['bc']['inlet']['p'][0]
        self.Dp = self._pOut - self._pIn
        L = -1
        convergence = [False, False]
        Qtol = 1E-10
        Q0old = self.Q0.copy()  # Old values, copied, otherwise we point to the same
        for region in self.regions:
            region.update()
        self.constants()  # Initialize the flow constants
        for i in range(50):
            # print('Starting initial condition iteration: ', i + 1)
            for r, region in enumerate(self.regions):
                tio = 0.5 * self._fluid['rho'] * (self._zetaIn / region.data['s'][0] ** 2 +
                                                  self._zetaOut / region.data['s'][L] ** 2)
                tcv = self.Wc(r, 1, region)[L] + self.Wv(r, 1, region)[L]

                self.Q0[r] = np.sqrt(-self.Dp / (tio + tcv))  # New Q0

                self.constants()  # Constants for the new Q0

                if (np.abs(self.Q0[r] - Q0old[r])) / np.abs(self.Q0[r]) < Qtol:
                    convergence[r] = True

            Q0old = self.Q0.copy()

            # print('Current Q0 is: ', self.Q0)
            # Check convergence of both regions
            if convergence[0] and convergence[1]:
                print("---> Initial Q0 has converged to: ", self.Q0, '. Iterations: ', i)
                break

        if not convergence[0] and not convergence[1]:
            raise ValueError("     ERROR: The fluid initial condition has not converged!")

    def update(self,  time, state):
        self._ti = time
        self.Q0 = state # Update the flow rate
        self.constants()
        # Update the regions and the operators
        for i, region in enumerate(self.regions):
            self.Q[i] = self.Q0[i] - region.data['dsi']
            self.v0[i] = self.Q0[i] / region.data['s'][0]

    def updateForces(self, time):
        # Update the boundary conditions
        if self._control['bc']['type'] == "variableInletPressure":
            pDict = self._control['bc']
            tpf = ((time - self._execution['time']['startTime']) /
                   (self._execution['time']['endTime'] - self._execution['time']['startTime']))
            self._pIn = pDict['inlet']['p'][0] + tpf * (pDict['inlet']['p'][-1] - pDict['inlet']['p'][0])

            self._pOut = pDict['outlet']['p']
            self.Dp = self._pOut - self._pIn

        # Force calculation for Calculix coupling
        # Integrate the pressure to obtain the force on every node
        F_half = 0.5 * si.cumtrapz(self.F(), self._mesh.x, initial=0.0)
        self.Forces[0:self._mesh.size, 1] = F_half
        self.Forces[self._mesh.size:self._mesh.size*2, 1] = F_half

    # Evaluate the RHS ofthe equation
    def rhs(self, time, state):
        rhs = np.zeros(self.dof)
        L = -1  # Index of the end of the channel
        self.update(None, state)
        # Assemble the rhs of each region
        deltaPx = np.zeros(1, dtype=object) # Pressure difference between the regions
        for i, region in enumerate(self.regions):
            s = region.data['s']
            dsi = region.data['dsi']

            t2 = 0.5 * self._fluid['rho'] * (
                      self._zetaIn * (self.Q[i][0] / s[0]) ** 2
                      + self._zetaOut * (self.Q[i][L] / s[L]) ** 2)

            Wcv = (self.Wc(i, self.Q[i]**2, region)
                  + self.Wv(i, self.Q[i]**2, region))

            t4 = - self.Wt(i, region.data['ddsi'], region)[L]  # Old value of acceleration

            rhs[i] = -(1 / self.Wt(i, 1, region)[L]) * (self.Dp + t2 + Wcv[L] + t4)

            # Update the acceleration (is this, the RHS)
            self.dQ0[i] = rhs[i]
            self.dQ[i] = self.dQ0[i] - region.data['ddsi']
            # Correct the pressure for the new acceleration
            self.px[i] = (self._pIn
                          - self._fluid['rho'] * 0.5 * self._zetaIn * (self.Q0[i] / region.data['s'][0]) ** 2
                          - self.Wt(i, self.dQ[i], region)
                          - Wcv)
            # print('exit pressure: ', self.px[i][-1])
            if region.type == 'top':
                sign = -1
            else:
                sign = 1
            deltaPx[0] += sign * self.px[i]

        # Update the pressure difference between the top and bottom
        self.deltaPx = deltaPx

        return rhs

    # Fluid pressure difference force
    def F(self):
        F = np.zeros(1, dtype=object)
        # for i, region in enumerate(self.regions):
        #     if region.type == 'top':
        #         sign = -1
        #     else:
        #         sign =  1
        #     F[0] += sign * self.px[i]
        F = -(self.px[0] - self.px[1])
        return F

    # ----- Flow Operators ----- #
    # Transient Operator
    def Wt(self, i, fx, region):
        size = region.data['s']
        Wt = self._fluid['rho'] * si.cumtrapz(fx / size, self._mesh.x, initial=0)
        return Wt

    # Convective Operator
    def Wc(self, i, fx, region):
        size = region.data['s']
        wc = si.cumtrapz((1/size) * np.gradient(fx / size, self._mesh.x, edge_order=2),
                          self._mesh.x, initial=0)
        Wc = self._fluid['rho'] * self.xix[i] * wc
        return Wc

    # Viscous Operator
    def Wv(self, i, fx, region):
        size = region.data['s']
        Wv = 0.25 * self.f0[i] * si.cumtrapz(fx / size**3, self._mesh.x, initial=0)
        return Wv

    # Write files
    def write(self):
        self.output[0].write(" ".join(map(str, self.Q0)) + '\n')
        self.output[1].write(" ".join(map(str, self.dQ0)) + '\n')
        # self.output[2].write(" ".join(map(str, self.px[0])) + '\n')
        # self.output[3].write(" ".join(map(str, self.px[1])) + '\n')
        self.output[2].write(" ".join(map(str, [self.px[0][0], self.px[1][0]])) + '\n')
        self.output[3].write(str(self._pIn) + '\n')
        self.output[4].write(" ".join(map(str, self.v0)) + '\n')
        self.output[5].write(str(self._ti) + '\n')
        # force0 = si.simps(-self.px[0] , self._mesh.x)
        # force1 = si.simps(self.px[1], self._mesh.x)
        # force = si.simps(-self.px[0] + self.px[1], self._mesh.x)
        # self.output[6].write(str(force) + '\n')

    # Calculate some constants
    def constants(self):
        for i, region in enumerate(self.regions):
            self.dRef[i] = region.data['s'][0]  # Reference channel size
            self.eRef[i] = self.dRef[i] / self.lRef
            self.vRef[i] = self.Q0[i] / self.dRef[i]
            self.Rd[i] = ReynoldsNumber(self, L=self.dRef[i], V=self.vRef[i])
            if self.Rd[i].value < 1000:  # Laminar
                if self.Rd[i].value <= 1: # Set the no flow values to zero to avoid indetermination
                    self.f0[i] = 0
                    self.xix[i] = 0
                    # self.eta = -0
                else: # Laminar values
                    self.f0[i] = 48.0 / self.Rd[i].value
                    self.xix[i] = 6.0 / 5.0
                    # self.eta = -self.f0 / Q0
            else:  # Turbulent
                self.f0[i] = 0.26 * self.Rd[i].value ** -0.24
                self.xix[i] = 1.0
                # self.eta = -(0.0624 * self.Rd.value ** -0.24) / Q0
