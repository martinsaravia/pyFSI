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
import scipy.integrate as si
import matplotlib.pyplot as plt

from pyFSI.vectors.eigen import eigenValue as eval, eigenVector as evec, eigenSystemVector as esys
import pyFSI.models.solidModels.calculixBeam as cx
from pyFSI.models.solidModels.solidBase import solidModel

class bernoulliEulerBeam(solidModel):
    def __repr__(self):
        return 'beamModel'

    def __init__(self, execution, control, mesh, time):
        super().__init__(execution, control, mesh, time)  # Call the base class

        # ----- Public attributes ----- #
        self.dof = control['solution']['modes']  # Modal DOFs
        self.m = None  # Mass vector
        self.k = None  # Stiffness vector
        self.c = None  # Damping vector
        self.K = np.zeros((self.dof, self.dof))  # Mass Matrix
        self.C = np.zeros((self.dof, self.dof))  # Stiffness Matrix
        self.M = np.zeros((self.dof, self.dof))  # Damping Matrix
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
        self.eigen = None  # Eigensystem
        # Dimensional analysis
        self.lRef = None  # Reference length
        self.tRef = None  # Reference thickness
        self.uRef = None  # Reference something
        self.vRef = None  # Reference velocity
        # Output variable mapping
        self.varMap["naturalFrequencies"] = "eigen.values"
        self.varMap["displacements"] = "y['mid']"
        self.varMap["velocities"] = "dy['mid']"
        self.varMap["accelerations"] = "ddy['mid']"
        self.varMap["numbers"] = "dimNumbers"

        # ----- Private attributes ----- #
        # Flags
        self._updated = False
        # Messages
        if self._debug:
            print("     WARNING: Beam initial conditions set to zero. ")

        # ----- Procedures ----- #
        self._initialize(control)

        self._stateMatrixSystem()  # Calculate the state matrix S

        self._stateVectorSystem()  # Calculate vectors k, m, c and F

        self.setInitialConditions()

    # Set the initial conditions
    def setInitialConditions(self):
        # Initial state conditions
        self.update(np.zeros(self.dof),
                    np.zeros(self.dof),
                    np.zeros(self.dof))

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

        self._updated = True

    # Build the dynamic vectors
    def _stateVectorSystem(self):
        self.m = self._m()  # Mass vector
        self.k = self._k()  # Stiffness vector
        self.c = self._c()  # Damping vector
        self.f = self._f()  # Load vector

    # Build the state matrices and vectors
    def _stateMatrixSystem(self):
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
        modes = self.eigen.vectors
        for i in range(self.dof):
            f[i] = si.simps(force * modes[i], self._mesh.x)
        F[self.dof:self.sof] = np.dot(self.Minv, f)
        return F

    def _m(self):  # Mass vector
        m = self._control['m'] * self.eigen.vectors
        return m

    def _c(self):  # Damping vector
        c = np.empty(len(self.eigen.values), dtype=object)
        for i, val in enumerate(c):
            c[i] = self._control['m'] * (self._control['solution']['damping'][0] *
                                         self.eigen.values[i] * self.eigen.vectors[i])
        return c

    def _k(self):  # Stiffness vector
        k = self._control['material']['E'] * self._control['I'] * self.eigen.d4
        return k

    def _f(self):  # body load
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

        elif pars['bc']['type'] == "simpleSupport":
            betaL = np.array([3.141592653589793, 6.283185307179586, 9.42477796076938, 12.566370614359172, 15.707963267948966, 18.84955592153876])
            beta = betaL / L
            vector = lambda x, b: np.sin(b * x)

        elif pars['bc']['type'] == "clampedSimple":
            betaL = np.array([3.9265975952148438, 7.068580627441406, 10.210182189941406, 13.351768493652344, 14.137166976928711, 16.49335479736328])
            beta = betaL / L
            vector = lambda x, b: (np.cos(b * x) - np.cosh(b * x))-((np.cos(b * L)-np.cosh(b * L))/(np.sin(b * L)-np.sinh(b * L))) * (np.sin(b * L)-np.sinh(b * L))

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
            if 'normalize' in self._control['solution']:
                normMethod = self._control['solution']['normalize']
                print("--> Choosing ", normMethod, "normalization for the beam...")
                if normMethod == 'no':
                    normMethod = False
            else:
                normMethod = False

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
        # Choose the boundary condition
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

        elif type == "simpleSupport":
            for i in range(es.size):
                es.setV(0, i, 0.0)
                es.setD2(0, i, 0.0)
                es.setV(-1, i, 0.0)
                es.setD2(-1, i, 0.0)
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

        elif type == "clampedSimple":
            for i in range(es.size):
                es.setV(0, i, 0.0)
                es.setD1(0, i, 0.0)
                es.setV(-1, i, 0.0)
                es.setD2(-1, i, 0.0)
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
        # self.output[0].write(" ".join(map(str, self.y['top'])) + '\n')
        # self.output[1].write(" ".join(map(str, self.y['mid'])) + '\n')
        # self.output[2].write(" ".join(map(str, self.y['bot'])) + '\n')
        # self.output[3].write(" ".join(map(str, self.dy['mid'])) + '\n')
        # self.output[4].write(" ".join(map(str, self.ddy['mid'])) + '\n')
        # self.output[5].write(" ".join(map(str, self.a)) + '\n')
        # self.output[6].write(" ".join(map(str, self.da)) + '\n')
        # self.output[7].write(" ".join(map(str, self.dda)) + '\n')
        pass

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





