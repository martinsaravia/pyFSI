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
from mesh.region.fsiRegion1D import fsiRegion1D
from models.flowModels.flowModel import flowModel
from properties.boundaryLayer import boundaryLayer
from vectors.eigen.eigenVector import eigenVector

class nlLeakageFlow2D(flowModel, ABC):
    def __repr__(self):
        return 'leakageFlow2DModel '

    def __init__(self, execution, control, mesh, boundary,  name='NN'):
        super().__init__(execution, control, mesh, boundary)

        # ----- Public attributes ----- #
        self.dof = 2  # Number of Q equations (1 per region)
        self.Q0 = np.zeros(self.dof)
        self.dQ0 = np.zeros(self.dof)
        self.Q = np.zeros(self.dof, dtype=object)
        self.dQ = np.zeros(self.dof, dtype=object)
        self.px = np.zeros(self.dof, dtype=object)
        # self.eRef = np.zeros(self.dof)
        # self.bLayer = np.zeros(self.dof)  # Boundary layer

        # ----- Private attributes ----- #
        self._th = 1  # Thickness of the model
        self._pIn = None  # Inlet pressure
        self._pOut = None  # Outlet pressure
        self._zetaIn = control['bc']['inlet']['zeta']  # Inlet loss factor
        self._zetaOut = control['bc']['outlet']['zeta']  # Outlet loss factor
        # self._f0 = np.zeros(self.dof)  # Viscous friction factor
        # self._xix = np.zeros(self.dof)  # Nonlinear profile factor
        # self._eta = np.zeros(self.dof)  # Derivative of f(Qx) at qx0
        self._gDof = None  # Size of the eigensystem
        # Convective operators
        self._Wc1 = np.zeros(self.dof, dtype=object)
        self._Wcdsi = np.zeros(self.dof, dtype=object)
        self._Wcdsi2 = np.zeros(self.dof, dtype=object)
        # Viscous operators
        self._Wv1 = np.zeros(self.dof, dtype=object)
        self._Wvdsi = np.zeros(self.dof, dtype=object)
        self._Wvdsi2 = np.zeros(self.dof, dtype=object)
        # Transient operators
        self._Wt1 = np.zeros(self.dof, dtype=object)
        self._Wtddsi = np.zeros(self.dof, dtype=object)
        # dimension of the associated eigensystem
        self._gDof = self.regions[0].eigen().size

        #Output
        self.output = []
        self.output.append(open(self._execution['paths']['fluidPath'] / 'Q0.out', 'w'))
        self.output.append(open(self._execution['paths']['fluidPath'] / 'dQ0.out', 'w'))

        # ----- Procedures ----- #
        # # Create the region objects
        # for i in control['regions']:
        #     region = fsiRegion1D(i, mesh, boundary)
        #     self.regions.append(region)
        # # dimension of the associated eigensystem
        # self._gDof = self.regions[0].eigen().size

        # Initial conditions
        # self.update(self.Q0, self.dQ0)

    # Flow rate equation initial condition
    def setInitialConditions(self):
        L = -1
        convergence = [False, False]
        Qtol = 1E-10
        print("Starting Q0 iteration with initial Q0 = ", self.Q0)
        Q0root = [0, 0]
        Q0pp = [0, 0]  # Parallel plates solution
        self.update(self.Q0, self.dQ0)
        for i in range(20):
            print('Starting initial condition iteration: ', i + 1)
            # Must be reset every iteration, otherwise self._Q0 is tied to it
            # from self._update
            Q0root = [0, 0]
            f0root = [0, 0] # Friction factor
            Q0pp = [0, 0] # Parallel plates solution
            for r, region in enumerate(self.regions):
                tio = 0.5 * self._fluid['rho'] * (self._zetaIn  / region.data['s'][0] ** 2 +
                                                  self._zetaOut / region.data['s'][L] ** 2)

                tcv = self._Wc1[r][L] + self._Wv1[r][L]

                tdp = self._pOut - self._pIn

                Q0root[r] = np.sqrt( -tdp / (tio + tcv ))

                f0root[r] = region.f0

                # print('r: ', r,
                #       'c: ', self._Wc1[r][L],
                #       'v: ', self._Wv1[r][L],
                #       'io: ', tio)

                if (np.abs(Q0root[r] - self.Q0[r])) / np.abs(Q0root[r]) < Qtol:
                    convergence[r] = True


            # # Print info
            # print("Q0 for iteration ", i, " is: ", Q0root)
            # print("f for iteration ", i, " is: ", f0root)



            # Update the regions and operators
            self.update(Q0root, self.dQ0)

            # Check convergence of both regions
            if convergence[0] and convergence[1]:
                print("---> Initial Q0 has converged to: ", Q0root)

                #  Inform the parallel plate solution
                for r, region in enumerate(self.regions):
                    Qpp = -(tdp * self._mesh.L) * region.data['s'][0] ** 3 / (
                                12 * self._fluid['mu'])
                    Q0pp[r] = Qpp
                print("The parallel plate solution is: ", Q0pp)

                break

    def update(self, Q0, dQ0):
        self.Q0 = Q0  # Update the flow rate
        self.dQ0 = dQ0  # Update the flow rate derivative
        # Update the boundary conditions
        if self._control['bc']['type'] == "variableInletPressure":
            tpfi = self._execution['time']['tpfi']
            self._pIn = (self._control['bc']['inlet']['p'][0] +
                         tpfi * (self._control['bc']['inlet']['p'][1] -
                                 self._control['bc']['inlet']['p'][0]))
            self._pOut = self._control['bc']['outlet']['p']

        # Update the regions and the operators
        for i, region in enumerate(self.regions):
            region.update(self, i)  # Update the region data
            self.Q[i] = self.Q0[i] - region.data['dsi']
            print('Qxt for region ', region.name, ' is: ', self.Q[i])
            self.dQ[i] = self.dQ0[i] - region.data['ddsi']
            # Convective operators
            self._Wc1[i] = self.Wc(1, region)
            self._Wcdsi[i] = self.Wc(region.data['dsi'], region)
            self._Wcdsi2[i] = self.Wc(region.data['dsi'] ** 2, region)
            # Viscous operators
            self._Wv1[i] = self.Wv(1, region)
            self._Wvdsi[i] = self.Wv(region.data['dsi'], region)
            self._Wvdsi2[i] = self.Wv(region.data['dsi'] ** 2, region)
            # Transient operators
            self._Wt1[i] = self.Wt(1, region)
            self._Wtddsi[i] = self.Wt(region.data['ddsi'], region)

            # Pressure equation
            self.px[i] = (self._pIn
                          - self._fluid['rho'] * 0.5 * self._zetaIn * (self.Q0[i] / region.data['s'][0]) ** 2
                          - self.Wt(self.dQ[i], region)
                          - self.Wv(self.Q[i]**2, region)
                          - self.Wc(self.Q[i]**2, region))

            # print('CHENNEL: ', i)
            # print('Q0: ', self.Q0[i])
            # print('s0: ', region.data['s'][0])
            # print('pin:', self._pIn)
            #
            # print('rhs_io: ', self._fluid['rho'] * 0.5 * self._zetaIn * (self.Q0[i] / region.data['s'][0]) ** 2)
            # print('Wt: ', self.Wt(self.dQ[i], region))
            # print('Wc: ', self.Wc(self.Q[i]**2, region))
            # print('Wv: ', self.Wv(self.Q[i]**2, region))

        print('PxTop: ', self.px[0])
        print('PxBot: ', self.px[1])

    # Evaluate the RHS ofthe equation
    def rhs(self):
        # Aliases
        rhs = np.zeros(2)
        L = -1  # Index of the end of the channel

        for i, region in enumerate(self.regions):
            s = region.data['s']
            dsi = region.data['dsi']

            # t1 = self._pOut - self._pIn
            # #
            # t2 =  ( 0.5 * self._fluid['rho'] * self._zetaIn * (1/ s[0])**2
            #       + 0.5 * self._fluid['rho'] * self._zetaOut * (1/ s[L])**2
            #       + self._Wc1[i][L]
            #       + self._Wv1[i][L]) * self.Q0[i]**2
            #
            # t3 = -(self._zetaOut * self._fluid['rho'] * dsi[L] / s[L]**2
            #       + 2 * self._Wcdsi[i][L]
            #       + 2 * self._Wvdsi[i][L]) * self.Q0[i]
            #
            # t4 = (0.5 * self._zetaOut * self._fluid['rho'] * (dsi[L] / s[L])**2
            #       - self._Wtddsi[i][L]
            #       + self._Wcdsi2[i][L]
            #       + self._Wvdsi2[i][L])
            #
            # rhs[i] = -(1 / self._Wt1[i][L]) * (t1 + t2 + t3 + t4)


            t1 = self._pOut - self._pIn

            t2 = 0.5 * self._fluid['rho'] * (
                      self._zetaIn * (self.Q[i][0] / s[0]) ** 2
                      + self._zetaOut * (self.Q[i][L] / s[L]) ** 2)

            t3 = self.Wc(self.Q[i]**2, region)[L] + self.Wv(self.Q[i]**2, region)[L]

            t4 = - self.Wt(region.data['ddsi'], region)[L]

            rhs[i] = -(1 / self._Wt1[i][L]) * (t1 + t2 + t3 + t4)

        return rhs

    # Added Mass Matrix (sum of the two regions)
    def M(self):
        M = np.zeros(self._gDof, dtype=object)
        L = -1
        # Add the added mass of each region
        for region in self.regions:  # Loop regions
            if region.type == 'top':
                sign = 1.0
            else:
                sign = -1.0
            vix = region.eigen().ix  # Indefinite integrals of the eigenvectors
            Wt1 = self.Wt(1, region)
            for i, ix in enumerate(vix):
                Wtg = self.Wt(ix, region)
                M[i] += -eigenVector(Wtg * (1 - Wt1 / Wt1[L]),
                                     self._mesh,
                                     normalize=False)

        return M

    # Added force
    def F(self):
        F = np.zeros(1, dtype=object)
        for i, region in enumerate(self.regions):
            if region.type == 'top':
                sign = 1.0
            else:
                sign = -1.0

            F[0] += -sign * self.px[i]

            # t1 = self._Wt1[i] * self.dQ0[i]

            # t1 = self.Wt(self.dQ[i], region)  # Este termino asi me implica sacar la matriz de masa agregada

            # t1mass = self.Wt1[i] * self.dQ0[i]
            #
            # t2 = self.Wc(self.Q[i]**2, region)  # Must be after calling self.dQ0
            #
            # t3 = self.Wv(self.Q[i]**2, region)
            #
            # t4 = self._fluid['rho'] * 0.5 * self._zetaIn * (self.Q0[i] / region.data['s'][0])**2
            #
            # # t4 = np.ones(self._mesh.size) * t4
            #
            # pxt = self._control['bc']['inlet']['p'][0] - (t1mass + t2 + t3 + t4)
            # print('pxt: ', pxt)
            #
            # print('Qxt: ', self.Q[i])
            # print('t1: ', t1, 't2: ', t2,'t3: ', t3,'t4: ', t4,)
            #
            # F[0] += sign * (t1 + t2 + t3 + t4)

        return F

    # ----- Flow Operators ----- #
    # Transient Operator
    def Wt(self, fx, region):
        size = region.data['s']
        Wt = self._fluid['rho'] * si.cumtrapz(fx / size, self._mesh.x, initial=0)
        return Wt

    # Convective Operator
    def Wc(self, fx, region):
        size = region.data['s']
        # Incorrect ? Is inconsistent with zero convection in parallel plates
        # c1 = fx / size**2
        # c2 = 0.5 * si.cumtrapz(fx * np.gradient(1 / size**2,
        #                                         self._mesh.x,
        #                                         edge_order=2),
        #                        self._mesh.x,
        #                        initial=0)
        # Wc = self._fluid['rho'] * region.xix * (c1 - c2)

        # Form without using the integral division rule
        cnew = si.cumtrapz( (1/size) * np.gradient(fx / size,
                                                self._mesh.x,
                                                edge_order=2),
                               self._mesh.x,
                               initial=0)
        Wc2 = self._fluid['rho'] * region.xix * cnew

        return Wc2

    # Viscous Operator
    def Wv(self, fx, region):
        size = region.data['s']
        Wv = 0.25 * region.f0 * si.cumtrapz(fx / size**3, self._mesh.x, initial=0)
        # Wv = 12 * self._fluid["mu"] * si.cumtrapz(fx / size ** 3, self._mesh.x, initial=0)
        return Wv


    def write(self):
        self.output[0].write(" ".join(map(str, self.Q0)) + '\n')
        self.output[1].write(" ".join(map(str, self.dQ0)) + '\n')



    


