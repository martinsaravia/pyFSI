from models.fluidModels.fluidModel import *
from cfd.dimensionlessNumbers.dimensionlessNumber import *
import scipy.integrate as si
import numpy as np
import sys


class leakageFlow2D(fluidModel):

    def __repr__(self):
        return 'leakageFlow2DModel '

    def __init__(self, control, mesh, boundary, dict, name='NN'):
        if control['debug'] == "y":
            print("----------------------" + self.__repr__() + " messages... ----------------------")
            print("--> Warning: Local Re is the same for all regions. Not as above Eq. 2.30")
            print("--> Warning: eta calculation must be checked")
            print("--> Warning: Check the he consideration of the beam thickness")
            print("--------------------------------------------------------------------------------")

        # References
        self._regions = dict['regions'] #Regions dict
        self._mesh = mesh
        self._bc = dict['bc']
        self._th = dict['thickness']

        # Flow rate (reads the key Qt that is added to the caseDict by the solver)
        self._qx0 = dict['bc']['inlet']['Qt']

        # Init the boundary related parameters
        hTop = boundary[dict['regions']['channelTop']['topBoundary']].ybot()[0]
        hBot = boundary[dict['regions']['channelBot']['botBoundary']].ytop()[0]
        self._area = dict['thickness'] * (hTop - hBot)
        #self._perimeter = 2 * dict['thickness'] + 2 *(hTop - hBot)
        #hydraulicDiameter = self._area / self._perimeter
        inletVelocity = self._qx0 / (dict['thickness'] * (hTop - hBot))

        # Init the base class
        self._Re = Reynolds((hTop - hBot),
                            inletVelocity,
                            dict['fluid']['type'],
                            dict['units'])
        super(leakageFlow2D, self).__init__(self._Re.fluid())

        # Create the region dictionary with useful things
        # All regiones must have at least/most one flexible boundary
        for key, val in self._regions.items():  # val is the boundary object
            val['name'] = key
            val['he'] = (boundary[val['topBoundary']].ybot() -
                         boundary[val['botBoundary']].ytop())
            val['he2'] = val['he']**2
            val['he3'] = val['he']**3
            val['he4'] = val['he']**4
            val['hed'] = np.gradient(val['he'], self._mesh.x(), edge_order=2)
            val['hex'] = si.cumtrapz(1.0/val['he'], self._mesh.x(), initial=0.0)
            val['heL'] = si.simps(1.0/val['he'], self._mesh.x())
            val['hexL'] = val['hex'] / val['heL']
            if boundary[val['topBoundary']].isFlexible(): # reference to the flexible boundary eigensystem
                val['eigen'] = boundary[val['topBoundary']].eigen()
            elif boundary[val['botBoundary']].isFlexible():
                val['eigen'] = boundary[val['botBoundary']].eigen()
            else:
                print("!!! WARNING: No flexible boundaries found in region" + rd['name'])

        # Nonlinear profile (xix), friction (ff) and inlet-outlet factors (eta)
        ReTol = 1E-10
        if self._Re.value() < 1:
            self._f0 = 48 * (1.0/(self._Re.value() + ReTol))
            self._xix = 6.0 / 5
            self._eta = 48.0 * self._nu * np.log(self._qx0 + ReTol)
        else:
            self._f0 = 0.26 * (self._qx0 / self._nu)**(-0.24)  # Friction factor
            self._xix = 1.0
            # self._eta = 0.26 * (self._nu**-0.24) * 1.31579 * self._qx0**(19.0/25)
            self._eta = -0.0624 / (self._nu * (self._qx0 / self._nu)**1.24)

        # Size of the system
        values = list(self._regions.values())
        self._size = values[0]['eigen'].dict()['size'] # Get the size from any of the regions

    def M(self):
        m = np.zeros(self._size, dtype=object)
        # Loop each region
        for r in self._regions.values():  # Loop for filling the mass matrix
            # Loop each mode
            for i in range(0, self._size):
                m[i] += self.integrate(r['eigen'].ix()[i] / r['he'], r)

        M = -self._rho * m

        return self._th * M

    def K(self):
        knl = np.zeros(self._size, dtype=object)  # Non Linear added stiffness
        kio = np.zeros(self._size, dtype=object)  # Input-Output added stiffness
        kvf = np.zeros(self._size, dtype=object)  # Viscous friction added stiffness

        # Loop each region
        for r in self._regions.values():  # Loop for filling the mass matrix
            # Loop each mode
            for i in range(0, self._size):
                gi = r['eigen'].v()[i]
                # Nonlinear profile stiffness
                tnl = 3 * gi * r['hed'] / r['he4'] - gi.d1() / r['he3']
                knl[i] += self._xix * self.integrate(tnl, r)
                # Viscous friction stiffness
                tvf = - gi / r['he4']
                kvf[i] += (3.0 * self._f0 / 4.0) * self.integrate(tvf, r)
                # Inlet-outlet stiffness
                tio0 = self._bc['inlet']['zeta'] * gi[0] / r['he3'][0]
                tio1 = self._bc['outlet']['zeta'] * gi[-1] / r['he3'][-1]
                kio[i] += (tio0 + tio1) * r['hexL'] - tio0

        K = self._rho * self._qx0**2 * (knl + kio + kvf)

        return self._th * K


    def C(self):
        cnl = np.zeros(self._size, dtype=object) # Non Linear added stiffness
        cio = np.zeros(self._size, dtype=object) # Input-Output added stiffness
        cvf = np.zeros(self._size, dtype=object) # Viscous friction added stiffness
        # Loop (and SUM!)each region
        for r in self._regions.values(): # Loop for filling the damping matrix
            x = self._mesh.x()
            # Loop each mode
            for i in range(0, self._size):
                gi = r['eigen'].v()[i]
                gix = r['eigen'].ix()[i]
                # Nonlinear profile damping
                tnl = gix * r['hed'] / r['he3'] - gi / r['he2']
                cnl[i] += 2.0 * self._xix * self.integrate(tnl, r)
                # Viscous friction damping
                tvf = - gix / r['he3']
                cvf[i] += (0.5 * self._f0 + 0.25 * self._qx0 * self._eta) * self.integrate(tvf, r)
                # Inlet-outlet damping
                cio[i] += (self._bc['outlet']['zeta'] / r['he2'][-1]) * r['eigen'].iL()[i] * r['hexL']

        C = self._rho * self._qx0 * (cnl + cio + cvf)

        return self._th * C

    def Tf(self):
        Tf = {}  # One vector T for each region
        for r in self._regions.values():
            tnl = np.zeros(1, dtype=object)
            tio = np.zeros(1, dtype=object)
            tvf = np.zeros(1, dtype=object)
            # Nonlinear term
            tnl = 2.0 * self._xix * self.integrate(r['hed'] / r['he3'], r)
            # Viscous Friction term
            tfv = -(0.5 * self._f0 + 0.25 * self._qx0 * self._eta) * self.integrate(1 / r['he3'], r)
            # Inlet-outlet Term
            t0 = self._bc['inlet']['zeta'] / r['he2'][0]
            t1 = self._bc['outlet']['zeta'] / r['he2'][-1]
            tio = (t0 + t1) * r['hexL'] - t0

            # Output as dictionary, one vector per region
            Tf[r['name']] = self._th * self._rho * self._qx0 * (tnl + tio + tvf)

        return Tf

    def Bq(self):
        Bq = {}
        for r in self._regions.values():  # Loop for filling the damping matrix
            Bq[r['name']] = np.zeros(1, dtype=object)
            for i in range(0, self._size):
                Bq[r['name']] += - self._th * si.simps(r['eigen'].ix()[i] / r['he'], self._mesh.x()) / r['heL']
        return Bq

    def Dq(self):
        Dq = {}
        for r in self._regions.values(): # Loop for filling the damping matrix
            dnl = np.zeros(1, dtype=object)
            dio = np.zeros(1, dtype=object)
            dvf = np.zeros(1, dtype=object)
            for i in range(0, self._size):
                # Nonlinear term
                dnl += 2 * self._xix * si.simps(r['eigen'].ix()[i] * r['hed'] / r['he3']
                    - r['eigen'].v()[i] / r['he2'], self._mesh.x())
                # Friction term
                dvf += -(0.5 * self._f0 + 0.25 * self._qx0 * self._eta) * si.simps(r['eigen'].ix()[i] / r['he3'], self._mesh.x())
                # Boundary Terms
                dio += -self._bc['outlet']['zeta'] / r['he2'][-1] * r['eigen'].iL()[i]

            Dq[r['name']] = self._th * self._qx0 / r['heL'] * (dnl + dio + dvf)

        return Dq

    def Eq(self):
        Eq = {}
        for r in self._regions.values(): # Loop for filling the damping matrix
            enl = np.zeros(1, dtype=object)
            eio = np.zeros(1, dtype=object)
            evf = np.zeros(1, dtype=object)
            for i in range(0, self._size):
                gi = r['eigen'].v()[i]
                gid = r['eigen'].d1()[i]
                # Nonlinear term
                enl += self._xix * si.simps(3 * gi * r['hed'] / r['he4'] - gid / r['he3'], self._mesh.x())
                # [r['name']]Friction Term
                evf += -(3.0 * self._f0 / 4) * si.simps(gi / r['he4'], self._mesh.x())
                # Boundary Term
                eio += -( self._bc['inlet']['zeta'] * gi[0] / r['he3'][0] + self._bc['outlet']['zeta'] * gi[-1] / r['he3'][-1] )

            Eq[r['name']] = self._th * self._qx0**2 / r['heL'] * (enl + eio + evf)

        return Eq

    def Gq(self):
        Gq = {}
        for r in self._regions.values(): # Loop for filling the damping matrix
            Gq[r['name']] = np.zeros(1, dtype=object)
            for i in range(0, self._size):
                gi = r['eigen'].v()[i]
                gid = r['eigen'].d1()[i]
                # Nonlinear term
                Gq[r['name']] += 2 * self._xix * si.simps(r['hed'] / r['he3'], self._mesh.x())
                # Friction Term
                Gq[r['name']] += -(0.5 * self._f0 + 0.25 * self._qx0 * self._eta) * si.simps(1 / r['he3'], self._mesh.x())
                # Boundary Term
                Gq[r['name']] += -(self._bc['inlet']['zeta'] / r['he2'][0] + self._bc['outlet']['zeta'] / r['he2'][-1] )

            Gq[r['name']] *= self._th * (self._qx0 / r['heL'])

        return Gq

    def integrate(self, term, region):
        # Routine for a common integration pattern in matrices from Tosi Appendix A
        x = self._mesh.x()
        return si.cumtrapz(term, x, initial=0.0) - si.simps(term, x) * region['hexL']

    def update(self, time):
        # Set the current q
        pass

    def xxqx0(self, region): # 2.59 Tosi
        if self._bc['inlet']['type']=='flowRate' or self._bc['outlet']['type']=='flowRate':
            print("--> ERROR: Flow rate can only be calculated with both pressure bcs.")
            sys.exit()
        pi = self._bc['inlet']['p']
        po = self._bc['outlet']['p']
        zetai = self._bc['inlet']['zeta']
        zetao = self._bc['outlet']['zeta']
        he0 = region['he'][0]
        heL = region['he'][-1]
        he = region['he']
        # Terms associated to nonlinear profile 2.23, input-output loss 2.45 and
        # viscous friction loss 2.29
        nlLoss = self._xix * 0.5 * (1/heL**2 - 1/he0**2)
        ioLoss = zetao / (2.0 * heL**2) + zetai / (2.0 * he0**2)
        vfLoss = 0.25 * self._f0 * si.simps(1 / he**3, self._mesh.x())

        self._qx0 = ( (pi - po) / (self._rho * (nlLoss + ioLoss + vfLoss)) )**0.5


    # Setters
    def setQ(self, val=0):
        self._qx0 = val


    # Getters
    def Q():
        return self._qx0

    def regions(self):
        return self._regions
        # for i, reg in enumerate(self._regions):
        #     mfi.append(self._rho * (integrate.simp(mesh.x(), ))
