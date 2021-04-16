# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: Leakage flow model of Tosi
#    I    #
# --------------------------------------------------------------------------- #
# Notes:
#   This class generates flow model based on Tosi's formulation
# Warnings:
#   Region data is handled through to the self.data dictionary
#   Thickness is set to 1
#   Local Re is the same for all regions. Not as above Eq. 2.30
#   eta calculation must be checked
#   Check the he consideration of the beam thickness
#   Only symmetric regions allowed (a single Q0)
# --------------------------------------------------------------------------- #
import sys
import numpy as np
import scipy.integrate as si
from pyFSI.mesh.region.fsiRegion1D import fsiRegion1D
from pyFSI.models.flowModels.flowBase import flowModel
from pyFSI.models.properties.boundaryLayer import boundaryLayer
from pyFSI.models.properties.dimensionlessNumbers import dimensionlessNumber
from pyFSI.fields.boundary import boundaryConditions


class leakageFlow2D(flowModel):
    def __repr__(self):
        return 'leakageFlow2DModel '

    def __init__(self, execution, control, mesh, boundary, time, name='NN'):

        super().__init__(execution, control, mesh, boundary, time)

        # ----- Public attribues ----- #
        self.Q0 = None  # Flow rate at (x,t) = (0, 0)
        self.v0 = None  # Flow speed
        self.eRef = None
        self.bLayer = None  # Boundary layer
        self.updated = False
        self.m = None  # Mass vector
        self.k = None  # Stiffness vector
        self.c = None  # Damping vector
        self.Tf = None
        self.Bq = None
        self.Dq = None
        self.Eq = None
        self.Gq = None
        # Output variable mapping
        self.varMap["flowRates"]  = "Q0"
        self.varMap["flowSpeeds"] = "v0"

        # ----- Private attributes ----- #
        self._th = 1
        self._f0 = None  # Viscous friction factor
        self._xix = None  # Nonlinear profile factor
        self._eta = None  # Derivative of f(Qx) at Q0
        self.dRef = None  # Channel inlet size reference length
        self.lRef = None  # Reference length
        self.vRef = None   # Reference velocity
        self.eRef = None   # Size quotient
        self._size = self.regions[0].eigen().size  # Size of the associated eigensystem

        # ----- Procedures ----- #
        # Initialize the boundary conditions
        # Flow Rate BC
        Q0BCDict = self._control['bc']['inlet']['Q0']
        self._Q0InBC = getattr(boundaryConditions, Q0BCDict['type'])(Q0BCDict)
        Q0BCDict = self._control['bc']['outlet']['Q0']
        self._Q0OutBC = getattr(boundaryConditions, Q0BCDict['type'])(Q0BCDict)
        # Loss factor BC
        zetaBCDict = self._control['bc']['inlet']['zeta']
        self._zetaInBC = getattr(boundaryConditions, zetaBCDict['type'])(zetaBCDict)
        zetaBCDict = self._control['bc']['outlet']['zeta']
        self._zetaOutBC = getattr(boundaryConditions, zetaBCDict['type'])(zetaBCDict)

        # Messages
        if self._debug:
            print("     WARNING: The region inlet size is assumed "
                  "to be the same for both regions. Check the flow model.")
            print("     WARNING: The inlet flow rate is the same for both regions.")

        self.update()

    def update(self):
        # Flow rate (reads the key Qt that is added to the casecontrol by the solver)
        self.Q0 = self._Q0InBC.getValue(self._time.value)

        # Update the geometry intermediate variables
        for region in self.regions:
            region.update()  # Update to read the initialized channel size
            # regionObj.data['name'] = region['name']
            region.data['he'] = region.data['s']   # Initial size of the region
            region.data['he0'] = region.data['he'][0]
            region.data['he2'] = region.data['he'] ** 2
            region.data['he3'] = region.data['he'] ** 3
            region.data['he4'] = region.data['he'] ** 4
            region.data['hed'] = np.gradient(region.data['he'], self._mesh.x, edge_order=2)
            region.data['hex'] = si.cumtrapz(1.0 / region.data['he'], self._mesh.x, initial=0.0)
            region.data['heL'] = si.simps(1.0 / region.data['he'], self._mesh.x)
            region.data['hexL'] = region.data['hex'] / region.data['heL']

        # Dimensionless numbers and boundary layer
        self.dRef = list(self.regions)[0].data['he0']
        self.lRef = self._mesh.L
        self.vRef = self.Q0 / self.dRef
        self.eRef = self.dRef / self.lRef
        self.bLayer = boundaryLayer(self, xPosition=self.lRef)
        self.v0 = self.Q0 / self.dRef
        super().calcNumbers()

        # Nonlinear profile (xix), viscous friction (f0) and derivative of f0 (eta)
        Rd = self.dimNumbers['Rd'].value
        if Rd < 1:  # Laminar
            self._f0 = 48.0 / Rd
            self._xix = 6.0 / 5.0
            self._eta = -self._f0 / self.Q0
        else:  # Turbulent
            self._f0 = 0.26 * Rd ** -0.24
            self._xix = 1.0
            self._eta = -(0.0624 * Rd**-0.24) / self.Q0
            # self._eta2 = -0.0624 * self._fluid['nu']**0.24 / self.Q0**1.24

        # Update the flow vectors
        self.m = self._m()
        self.k = self._k()
        self.c = self._c()
        self.Tf = self._Tf()
        self.Bq = self._Bq()
        self.Dq = self._Dq()
        self.Eq = self._Eq()
        self.Gq = self._Gq()

        self.updated = True

    # Mass vector
    def _m(self):
        m = np.zeros(self._size, dtype=object)
        # Loop each region
        for region in self.regions:  # Loop for filling the mass matrix
            r = region.data
            # Loop each mode
            for i in range(0, self._size):
                m[i] += self.integrate(region.eigen().ix[i] / r['he'], region)

        M = -self._fluid['rho'] * m

        return self._th * M

    # Stiffness vector
    def _k(self):
        knl = np.zeros(self._size, dtype=object)  # Non Linear added stiffness
        kio = np.zeros(self._size, dtype=object)  # Input-Output added stiffness
        kvf = np.zeros(self._size, dtype=object)  # Viscous friction added stiffness

        # Loop each region
        for region in self.regions:  # Loop for filling the mass matrix
            # Loop each mode
            for i in range(0, self._size):
                gi = region.eigen().vectors[i]
                r = region.data
                # Nonlinear profile stiffness
                tnl = 3 * gi * r['hed'] / r['he4'] - gi.d1 / r['he3']
                knl[i] += self._xix * self.integrate(tnl, region)
                # Viscous friction stiffness
                tvf = - gi / r['he4']
                kvf[i] += (3.0 * self._f0 / 4.0) * self.integrate(tvf, region)
                # Inlet-outlet stiffness
                tio0 = self._zetaInBC.getValue(self._time.value) * gi[0] / r['he3'][0]
                tio1 = self._zetaOutBC.getValue(self._time.value) * gi[-1] / r['he3'][-1]
                kio[i] += (tio0 + tio1) * r['hexL'] - tio0

        K = self._fluid['rho'] * self.Q0**2 * (knl + kio + kvf)

        return self._th * K

    # Damping vector
    def _c(self):
        cnl = np.zeros(self._size, dtype=object)  # Non Linear added stiffness
        cio = np.zeros(self._size, dtype=object)  # Input-Output added stiffness
        cvf = np.zeros(self._size, dtype=object)  # Viscous friction added stiffness
        # Loop (and SUM!)each region
        for region in self.regions:  # Loop for filling the damping matrix
            x = self._mesh.x
            # Loop each mode
            for i in range(self._size):
                gi = region.eigen().vectors[i]
                r = region.data
                gix = region.eigen().ix[i]
                # Nonlinear profile damping
                tnl = gix * r['hed'] / r['he3'] - gi / r['he2']
                cnl[i] += 2.0 * self._xix * self.integrate(tnl, region)
                # Viscous friction damping
                tvf = - gix / r['he3']
                cvf[i] += ((0.5 * self._f0 + 0.25 * self.Q0 * self._eta) * 
                           self.integrate(tvf, region))
                # Inlet-outlet damping
                cio[i] += ((self._zetaOutBC.getValue(self._time.value) / r['he2'][-1]) * 
                           region.eigen().iL[i] * r['hexL'])

        C = self._fluid['rho'] * self.Q0 * (cnl + cio + cvf)

        return self._th * C

    # Forcing vector
    def _Tf(self):
        Tf = {}  # One vector T for each region
        for region in self.regions:
            r = region.data
            tnl = np.zeros(1, dtype=object)
            tio = np.zeros(1, dtype=object)
            tvf = np.zeros(1, dtype=object)
            # Nonlinear term
            tnl = 2.0 * self._xix * self.integrate(r['hed'] / r['he3'], region)
            # Viscous Friction term
            tfv = -((0.5 * self._f0 + 0.25 * self.Q0 * self._eta) * 
                    self.integrate(1 / r['he3'], region))
            # Inlet-outlet Term
            t0 = self._zetaInBC.getValue(self._time.value) / r['he2'][0]
            t1 = self._zetaOutBC.getValue(self._time.value) / r['he2'][-1]
            tio = (t0 + t1) * r['hexL'] - t0
            # Output as dictionary, one vector per region
            Tf[region.name] = self._th * self._fluid['rho'] * self.Q0 * (tnl + tio + tvf)

        return Tf

    def _Bq(self):
        Bq = {}
        for region in self.regions:  # Loop for filling the damping matrix
            r = region.data
            Bq[region.name] = np.zeros(1, dtype=object)
            for i in range(self._size):
                Bq[region.name] += - (self._th * si.simps(region.eigen().ix[i] /
                                    r['he'], self._mesh.x) / r['heL'])
        return Bq

    def _Dq(self):
        Dq = {}
        for region in self.regions: # Loop for filling the damping matrix
            r = region.data
            dnl = np.zeros(1, dtype=object)
            dio = np.zeros(1, dtype=object)
            dvf = np.zeros(1, dtype=object)
            for i in range(self._size):
                # Nonlinear term
                dnl += (2 * self._xix * si.simps(region.eigen().ix[i] * r['hed'] / r['he3']
                        - region.eigen().vectors[i] / r['he2'], self._mesh.x))
                # Friction term
                dvf += -((0.5 * self._f0 + 0.25 * self.Q0 * self._eta) * 
                         si.simps(region.eigen().ix[i] / r['he3'], self._mesh.x))
                # Boundary Terms
                dio += -self._zetaOutBC.getValue(self._time.value) / r['he2'][-1] * region.eigen().iL[i]

            Dq[region.name] = self._th * self.Q0 / r['heL'] * (dnl + dio + dvf)

        return Dq

    def _Eq(self):
        Eq = {}
        for region in self.regions:  # Loop for filling the damping matrix
            r = region.data
            enl = np.zeros(1, dtype=object)
            eio = np.zeros(1, dtype=object)
            evf = np.zeros(1, dtype=object)
            for i in range(0, self._size):
                gi = region.eigen().vectors[i]
                gid = region.eigen().d1[i]
                # Nonlinear term
                enl += (self._xix * si.simps(3 * gi * r['hed'] / r['he4'] 
                        - gid / r['he3'], self._mesh.x))
                # [region.name]Friction Term
                evf += -(3.0 * self._f0 / 4) * si.simps(gi / r['he4'], self._mesh.x)
                # Boundary Term
                eio += -((self._zetaInBC.getValue(self._time.value) * gi[0] / r['he3'][0] 
                          + self._zetaOutBC.getValue(self._time.value) * gi[-1] / r['he3'][-1]))

            Eq[region.name] = self._th * self.Q0**2 / r['heL'] * (enl + eio + evf)

        return Eq

    def _Gq(self):
        Gq = {}
        for region in self.regions:  # Loop for filling the damping matrix
            r = region.data
            Gq[region.name] = np.zeros(1, dtype=object)
            for i in range(self._size):
                gi = region.eigen().vectors[i]
                gid = region.eigen().d1[i]
                # Nonlinear term
                Gq[region.name] += 2 * self._xix * si.simps(r['hed'] / r['he3'],
                                                            self._mesh.x)
                # Friction Term
                Gq[region.name] += -((0.5 * self._f0 + 0.25 * self.Q0 * self._eta) *
                                   si.simps(1 / r['he3'], self._mesh.x))
                # Boundary Term
                Gq[region.name] += -((self._zetaInBC.getValue(self._time.value) / r['he2'][0]
                                     + self._zetaOutBC.getValue(self._time.value) / r['he2'][-1]))

            Gq[region.name] *= self._th * (self.Q0 / r['heL'])

        return Gq

    def integrate(self, term, region):
        # Routine for a common inte
        # gration pattern in matrices from Tosi Appendix A
        x = self._mesh.x
        return (si.cumtrapz(term, x, initial=0.0) -
                si.simps(term, x) * region.data['hexL'])




