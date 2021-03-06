# --------------------------------------------------------------------------- #
#    p    #     version: 0.0
#    y    #     date: 20/03/2020
#    F    #     author: Martin Saravia
#    S    #     description: BE Beam class
#    I    #     return: solidModel object
# --------------------------------------------------------------------------- #
# Notes:
#   The class  works either as an analytical model or with a calculix model
#   The model type is selected from the key "method" in the solid-solution dict
# --------------------------------------------------------------------------- #
import numpy as np
import scipy.optimize as opt
from scipy import interpolate
from vectors.eigen import eigenValue as eval, eigenVector as evec, eigenSystemVector as esys
import models.solidModels.calculixBeam as cx
from models.solidModels.solidModel import solidModel

class bernoulliEulerBeam(solidModel):
    def __repr__(self):
        return 'beamModel'

    def __init__(self, execution, control, mesh, name='beam'):
        if execution['debug'] == "y":
            print("----------------------" + self.__repr__() + " messages... ----------------------")
            print("--> Warning: modal damping is equal for every mode...")
            print("--> Warning: reference displacments fixed to L/10...")
            print("----------------------------------------------------------------------------")

        super().__init__(execution, control, mesh, name)

        # Kinematics
        self._y = np.zeros(len(mesh.x()))  # Initial y set to zero
        self.xi = mesh.x()[0]  # Initial position
        self.xf = mesh.x()[-1]  # Final position
        self.L = self.xf - self.xi

        # Add calculated parameters to the input controlionary
        control['L'] = self.L  # Add the length to the parameters dictionary
        control['I'] = (control['section']['b'] * control['section']['h'] ** 3
                          / 12.0)
        control['A'] = control['section']['b'] * control['section']['h']
        control['m'] = control['A'] * control['material']['rho']

        # Create the eigensystem object
        if control['solution']['type'] == 'modal':
            if control['solution']['method'] == 'analytic':
                self._ES = self._analyticEigenSystem(control['solution']['modes'])
            elif control['solution']['method'] == 'calculix':
                self._ES = self._calculixEigenSystem(control['solution']['modes'])
            else:
                print("--> ERROR: Beam solution method " +
                      control['solution']['method']
                      + " not kwnown !")

            self._setBC(control['bc']['type'], self._ES)  # Impose the boundary condition

            self._ES.calculate()  # Calculate integrals

        else:
            print("--> ERROR: Beam solution type " + control['solution']['type'] + " not kwnown !")

        # Fill the reference parameters
        self.lRef = self.L
        self.uRef = self.L / 10.0
        # We take the largest natural frequency to calculate the reference velocity
        # Scaling with the sound speed in the solid would be not very realistic
        self.vRef = self.freqs()[-1] * self.uRef

        # Calculate the dimensionless numbers
        self.calcNumbers()

    # Mass vector
    def M(self):
        return self._control['m'] * self._ES.v()

    # Damping vector
    def C(self):
        damp = np.empty(len(self._ES.values()), dtype=object)
        for i, val in enumerate(damp):
            damp[i] = self._control['solution']['damping'] * self._ES.values()[i] * self._ES.v()[i]
        return self._control['m'] * damp

    # Stiffness vector
    def K(self):
        return self._control['material']['E'] * self._control['I'] * self._ES.d4()

    # Private functions
    def _analyticEigenSystem(self, nmodes):
        # Aliases
        mesh = self._mesh
        pars = self._control
        L = self.L

        # Select the type of bc
        if pars['bc']['type'] == "clampedFree":
            betaL = np.array(
                [1.875104068711961, 4.694091132974175, 7.854757438237612,
                 10.995540734875467, 14.13716839104647,17.278759532088234]
            )
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
            cc = ((pars['material']['E'] * pars["I"]) / (pars['material']['rho'] * pars["A"])) ** 0.5
            values.append(eval.eigenValue(beta[i] ** 2 * cc))
            # Eigenvectors
            vectors.append(evec.eigenVector(vector(mesh.x(), beta[i]), mesh, info="Beam eigenvector " + str(i + 1),
                                            normalize=False))

        # Create the eigensystem
        eigenSystem = esys.eigenSystemVector(values, vectors)

        return eigenSystem

    def _calculixEigenSystem(self, nmodes):

        ccxModel = cx.calculixBeam(self._control, self._mesh, self._control)

        return ccxModel.eigenSystem()

    def _setBC(self, type, eigenSystem):
        # Impose bcs
        es = self._ES
        mesh = self._mesh
        if type == "clampedFree":
            # Create the interpolate object
            for i in range(0, es.dict()['size']):
                es.setV(0, i, 0.0)
                es.setD1(0, i, 0.0)
                es.setD2(-1, i, 0.0)
                es.setD3(-1, i, 0.0)
                # Correct the third derivative
                itpobj = interpolate.interp1d(mesh.x()[3:], es.d3()[i][3:], fill_value='extrapolate')
                es.setD3(range(0, 3), i, itpobj(mesh.x()[0:3]))
                # Correct the fourth derivative
                itpobj = interpolate.interp1d(mesh.x()[4:-4], es.d4()[i][4:-4], fill_value='extrapolate')
                es.setD4(range(0, 4), i, itpobj(mesh.x()[0:4]))
                es.setD4(range(-4, 0), i, itpobj(mesh.x()[-4:]))

    def __eigenSolve(self, modes=50, freqs=10):
        betamin = 2
        betamax = freqs
        betatol = 1
        freqv = []
        # Find the betas and the modes
        # if self.bc == "clampedFree":
        #     charEqn = lambda beta : np.cosh(beta * self.L) * np.cos(beta * self.L) + 1.0
        mode = 0
        freq = 0
        while mode <= modes:
            # Characteristic equations and its derivative
            if self._control['bc']['type'] == "clampedFree":
                charEqn = lambda beta: np.cosh(beta * self.L) * np.cos(beta * self.L) + 1.0
                charEqnPrime = lambda beta: self.L * np.sinh(beta * self.L) * np.cos(beta * self.L) - self.L * np.sin(
                    beta * self.L) * np.cosh(beta * self.L)

            eval, conv = opt.newton(charEqn, betamin, fprime=charEqnPrime, full_output=True, maxiter=500)

            if conv.converged:
                betamin += betatol
                if abs(abs(eval) - abs(freq)) > betatol:
                    freq = eval
                    freqv.append(freq)
                    mode += 1
                    print("Mode " + str(mode) + " found at " + str(eval * self.L))
        # Filter and sort
        freqv = [item * self.L for item in freqv if item >= 0]
        freqv.sort()

    def charEqn(self, beta):
        f = None
        if self._control['bc']['type'] == "clampedFree":
            f = np.cosh(beta * self.L) * np.cos(beta * self.L) + 1.0
        return f

    # Getters
    def y(self):
        return self._y

    def ytop(self):
        return self._y + 0.5 * self._control['section']['h']

    def ybot(self):
        return self._y - 0.5 * self._control['section']['h']

    def eigen(self):
        return self._ES

    def mass(self):
        return self._control['A'] * self._control['material']['rho']

    def freqs(self):
        return self._ES.values()
