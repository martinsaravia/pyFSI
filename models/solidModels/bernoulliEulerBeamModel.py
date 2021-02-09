# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: BE Beam class
#    I    #     return: solidModel object
# --------------------------------------------------------------------------- #
# Notes:
#   The class  works either as an analytical model or with a Calculix model
#   The model type is selected from the key "method" in the solid-solution dict
#   Supports reconstruction of displacements
# Warnings:
#   Modal damping is equal for every mode
#   Reference displacments fixed to L/10
#   Initial position of the beam is set to zero
# Optimize:
#   Delete info from dy and ddy
# --------------------------------------------------------------------------- #

import numpy as np
import scipy.optimize as opt
from scipy import interpolate
from vectors.eigen import eigenValue as eval, eigenVector as evec, eigenSystemVector as esys
import models.solidModels.calculixBeam as cx
from models.solidModels.solidBase import solidModel
import scipy.integrate as si
import matplotlib.pyplot as plt

class bernoulliEulerBeam(solidModel):
    def __repr__(self):
        return 'beamModel'

    def __init__(self, execution, control, mesh, name='beam'):
        super().__init__(execution, control, mesh, name)  # Call the base class

        # ----- Public attributes ----- #
        self.dof = control['solution']['modes']  # Modal DOFs
        self.K = np.zeros((self.dof, self.dof))
        self.C = np.zeros((self.dof, self.dof))
        self.M = np.zeros((self.dof, self.dof))
        self.Minv = np.zeros((self.dof, self.dof))
        self.sof = 2 * self.dof
        self.S = np.zeros((self.sof, self.sof)) # State Matrix
        self.F = np.zeros(self.sof)  # Modal load
        self.state = np.zeros(self.sof)
        self.a = None  # Generalized displacement
        self.da = None  # Generalized velocity
        self.dda = None  # Generalized acceleration
        self.xi = None  # Initial position
        self.xf = None  # Final position
        self.y = {}   # Current position of top, bot and mid surfaces
        self.dy = {}  # Current velocity of top, bot and mid surfaces
        self.ddy = {}  # Current accelerations of top, bot and mid surfaces
        # Eigensystem
        self.eigen = None
        # Dimensional analysis
        self.lRef = None  # Reference length
        self.tRef = None  # Reference thickness
        self.uRef = None  # Reference something
        self.vRef = None  # Reference velocity

        # ----- Private attributes ----- #
        self._mesh = mesh

        # ----- Procedures ----- #
        self._initialize(control)

        self._stateSystem()

        self.setInitialConditions()


    # Set the initial conditions
    def setInitialConditions(self):
        # Initial state conditions
        self.update(np.zeros(self.dof),
                    np.zeros(self.dof),
                    np.zeros(self.dof))

    # ----- Public methods ----- #
    # Update the geometry
    def update(self, a, da, dda):
        self.a = a  # State variables
        self.da = da
        self.dda = dda
        # Update the displacements
        self.y['mid'] = self.eigen.reconstruct(self.a)
        self.y['top'] = self.y['mid'] + self._control['hm']
        self.y['bot'] = self.y['mid'] - self._control['hm']
        # Update the velocities
        self.dy['mid'] = self.eigen.reconstruct(self.da)
        # Update the accelerations
        self.ddy['mid'] = self.eigen.reconstruct(self.dda)

    # Build the state matrix
    def _stateSystem(self):
        # State matrix assumming mass normalized eigenvectors
        dof = self.dof
        sof = self.sof
        I = np.identity(dof)

        # Linear mass
        self.M = I
        self.Minv = self.M

        # Assemble the linear stiffness for self-adjoint operator
        for i in range(dof):
            self.K[i, i] = self.eigen.values[i] ** 2

        # Assemble the damping matrix
        if dof > 1:
            s1 =  self._control['solution']['damping'][0]
            s2 = self._control['solution']['damping'][1]
            w1 = self.eigen.values[0]
            w2 = self.eigen.values[1]
            alpha = 2 * ( s2 * w1**2 * w2 - s1 * w1 * w2**2) / (w1**2 - w2**2)
            beta = 2 * ( s1 * w1 - s2 * w2) / (w1**2 - w2**2)
            for i in range(dof):
                self.C[i, i] = alpha * self.M[i, i] + beta * self.K[i, i]
        else:
            self.C[0, 0] = 2 * self._control['solution']['damping'][0] * self.eigen.values[0]


        # Assemble the modal force
        f = np.zeros(dof)
        for i in range(dof):
            f[i] = self._control['w'] * si.simps(self.eigen.vectors[i], self._mesh.x)  # gravity
        self.F[dof:sof] = np.dot(self.Minv, f)

        # State Matrix
        self.S[0:dof, dof:sof] = I
        self.S[dof:sof, 0:dof] = -np.dot(self.Minv, self.K)
        self.S[dof:sof, dof:sof] = -np.dot(self.Minv, self.C)


    # Add a fluid force
    def addedStateModalForce(self, force):
        F = np.zeros(self.sof)
        f = np.zeros(self.dof)
        mode = self.eigen.vectors
        for i in range(self.dof):
            # plt.plot(force)
            # plt.figure()
            # plt.plot(force  * mode[i])
            # plt.show()
            f[i] = si.simps(force * mode[i], self._mesh.x)
        F[self.dof:self.sof] = np.dot(self.Minv, f)

        return F

    def m(self):  # Mass vector
        m = self._control['m'] * self.eigen.vectors
        return m

    def c(self):  # Damping vector
        c = np.empty(len(self.eigen.values), dtype=object)
        for i, val in enumerate(c):
            c[i] = self._control['m'] * (self._control['solution']['damping'][0] *
                                         self.eigen.values[i] * self.eigen.vectors[i])
        return c

    def k(self):  # Stiffness vector
        k = self._control['material']['E'] * self._control['I'] * self.eigen.d4
        return k

    def f(self):  # body load
        # F = self._control['w'] * np.ones(len(self.eigen.vectors[0]))
        f = np.zeros(1, dtype=object)
        f[0] = self._control['w'] * np.ones(len(self.eigen.vectors[0]))
        return f

    # ----- Private methods ----- #
    def _initialize(self, control):
        # Add calculated parameters to the control dictionary
        control['L'] = self._mesh.L  # Add the length to the parameters dictionary
        control['I'] = (control['section']['b'] * control['section']['h'] ** 3 / 12)
        control['A'] = control['section']['b'] * control['section']['h']
        control['m'] = control['A'] * control['material']['rho']  # mass per unit
        if 'g' in control["section"]:
            control['w'] = control['m'] * control['section']['g']  # weight per unit lenght
        else:
            control['w'] = 0.0
        control['hm'] = 0.5 * self._control['section']['h']  # Mid-distance

        # Create the eigensystem object
        if control['solution']['type'] == 'modal':
            if control['solution']['method'] == 'analytic':
                self.eigen = self._analyticEigenSystem(control['solution']['modes'])
            elif control['solution']['method'] == 'calculix':
                self.eigen = self._calculixEigenSystem(control['solution']['modes'])
            else:
                print("--> ERROR: Beam solution method " +
                      control['solution']['method'] + " not kwnown !")
            # Impose the boundary condition on the mode
            self._setBC(control['bc']['type'],
                        self.eigen)
            # Correct boundary conditions in the derivatives
            self.eigen.correct()

        else:
            print("--> ERROR: Beam solution type " + control['solution'][
                'type'] + " not kwnown !")


        # Fill the reference parameters
        self.lRef = control['L']
        self.tRef = control['section']['h']
        self.uRef = control['L'] / 10.0
        # We take the largest natural frequency to calculate the reference velocity
        # Scaling with the sound speed in the solid would be not very realistic
        self.vRef = self.freqs()[-1] * self.uRef

        # Calculate the dimensionless numbers
        self.calcNumbers()

        # Initialize the output files
        self.output.append(open(self._execution['paths']['solidPath'] / 'yTop.out', 'w'))
        self.output.append(open(self._execution['paths']['solidPath'] / 'yMid.out', 'w'))
        self.output.append(open(self._execution['paths']['solidPath'] / 'yBot.out', 'w'))
        self.output.append(open(self._execution['paths']['solidPath'] / 'dyMid.out', 'w'))
        self.output.append(open(self._execution['paths']['solidPath'] / 'ddyMid.out', 'w'))
        self.output.append(open(self._execution['paths']['solidPath'] / 'a.out', 'w'))
        self.output.append(open(self._execution['paths']['solidPath'] / 'da.out', 'w'))
        self.output.append(open(self._execution['paths']['solidPath'] / 'dda.out', 'w'))

        # Save the mesh
        self.output.append(open(self._execution['paths']['solidPath'] / 'x.out', 'w'))
        self.output[8].write(" ".join(map(str, self._mesh.x)) + '\n')
        self.output[8].close()

        # Initial conditions
        # y0 = np.zeros(len(self._mesh.x))  # Initial position
        # y0[:] = control['bc']['y0']
        # dy0 = np.zeros(len(self._mesh.x))  # Initial velocity
        # dy0[:] = control['bc']['dy0']
        # ddy0 = np.zeros(len(self._mesh.x))  # Zero initial acceleration

    def _analyticEigenSystem(self, nmodes):
        # Aliases
        mesh = self._mesh
        pars = self._control
        L = self._control['L']

        # Select the type of bc
        if pars['bc']['type'] == "clampedFree":
            betaL = np.array(
                [1.875104068711961, 4.694091132974175, 7.854757438237612,
                 10.995540734875467, 14.13716839104647,17.278759532088234])
            beta = betaL / L
            # The eigenvector function
            # vector = lambda x, b: np.cosh(b*x) - np.cos(b*x) + ((np.cos(b*L)+np.cosh(b*L))/(np.sin(b*L)+np.sinh(b*L)))*(np.sin(b*x)-np.sinh(b*x))
            vector = lambda x, b: np.cosh(b * x) - np.cos(b * x) - (
                        (np.sinh(b * L) - np.sin(b * L)) / (np.cosh(b * L) + np.cos(b * L))) * (
                         np.sinh(b * x) - np.sin(b * x))

        # Fill the eigensystem
        values = []
        vectors = []
        for i in range(0, nmodes):
            # Natural frequencies (not the eigenvalue beta)
            cc = ((pars['material']['E'] * pars["I"])
                  / (pars['material']['rho'] * pars["A"])) ** 0.5
            values.append(eval.eigenValue(beta[i] ** 2 * cc))
            # Eigenvectors

            # Choose the normalization method, No is for stiffness based on the
            # fourth derivative of the mode, mass is for k based on the eigenvalue
            if self._control['solution']['normalize'] == "No":
                normMethod = None
            elif self._control['solution']['normalize'] == "mass":
                print("--> Choosing mass normalization for the beam...")
                normMethod = "mass"
            vectors.append(evec.eigenVector(vector(mesh.x, beta[i]),
                                            mesh,
                                            info="Beam eigenvector " + str(i + 1),
                                            normalize=normMethod,
                                            mass=self._control['m']))

        # Create the eigensystem
        eigenSystem = esys.eigenSystemVector(values, vectors)
        return eigenSystem

    # Create an eigensystem based on a Calculix solution
    def _calculixEigenSystem(self, nmodes):
        ccxModel = cx.calculixBeam(self._control, self._mesh, self._control)
        return ccxModel.eigenSystem()

    # Set the boundary conditions on the eigenvectors
    def _setBC(self, type, eigenSystem):
        mesh = self._mesh
        es = self.eigen
        if type == "clampedFree":
            # Create the interpolate object
            for i in range(es.size):
                es.setV(0, i, 0.0)
                es.setD1(0, i, 0.0)
                es.setD2(-1, i, 0.0)
                es.setD3(-1, i, 0.0)
                # Correct the third derivative
                itpobj = interpolate.interp1d(mesh.x[3:],
                                              es.d3[i][3:],
                                              fill_value='extrapolate')
                es.setD3(range(0, 3), i, itpobj(mesh.x[0:3]))
                # Correct the fourth derivative
                itpobj = interpolate.interp1d(mesh.x[4:-4],
                                              es.d4[i][4:-4],
                                              fill_value='extrapolate')
                es.setD4(range(0, 4), i, itpobj(mesh.x[0:4]))
                es.setD4(range(-4, 0), i, itpobj(mesh.x[-4:]))

    def write(self):
        self.output[0].write(" ".join(map(str, self.y['top'])) + '\n')
        self.output[1].write(" ".join(map(str, self.y['mid'])) + '\n')
        self.output[2].write(" ".join(map(str, self.y['bot'])) + '\n')
        self.output[3].write(" ".join(map(str, self.dy['mid'])) + '\n')
        self.output[4].write(" ".join(map(str, self.ddy['mid'])) + '\n')
        self.output[5].write(" ".join(map(str, self.a)) + '\n')
        self.output[6].write(" ".join(map(str, self.da)) + '\n')
        self.output[7].write(" ".join(map(str, self.dda)) + '\n')

    # Getters
    def freqs(self):
        return self.eigen.values

    def mesh(self):
        return self._mesh

    def control(self):
        return self._control

    # Messages
    def _messages(self, n):
        pass





