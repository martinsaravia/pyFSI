import scipy.integrate as si
import numpy as np

from pyFSI.solvers.solverBase import solverBase


# Solver for transient FSI simulations
class transient(solverBase):
    def __init__(self, fsi):
        super().__init__(fsi)

    def solve(self):
        time = self._execution['time']
        initialConditions = self.fsi.state
        tspan = np.array([time['startTime'], time['startTime'] + time['deltaT']])
        dt = time['deltaT']

        # Choose integration scheme
        while tspan[1] <= time['endTime']:
            print("Solving for time: ", tspan[1])
            getattr(self, self.control['coupling'])(tspan)  # Call the implicit or the explicit step
            self.write(tspan[1])
            tspan += dt
        self.fsi.finish()

    def implicit(self, tspan):
        for n in range(1, self.control['subcycles'] + 1):
            pass

    def explicit(self, tspan):
        fsi = self.fsi

        # 1) Solve the flow
        # print("---> Solving the flow")
        fSol = si.solve_ivp(fsi.flow().rhs,
                            tspan,
                            fsi.fState,
                            method=self.control['integrator'][0],
                            atol=self.control['atol'],
                            rtol=self.control['rtol'])
        fsi.update('flow', tspan[1], fSol.y[:, -1])

        # 2) Update the boundary condition
        # print("---> Updating the flow force over the solid")
        fsi.update('forces', tspan[1], None)

        # 3) Solve the solid
        # print("---> Solving the solid")
        sSol = si.solve_ivp(fsi.rhsSolid,
                            tspan,
                            fsi.sState,
                            method=self.control['integrator'][1],
                            atol=self.control['atol'],
                            rtol=self.control['rtol'])
        fsi.update('solid', tspan[1], sSol.y[:, -1])

        # 4) Update the region geometry
        # print("---> Updating the region position, velocities and accelerations")
        fsi.update('regions', None, None)

    def write(self, t):
        # print("---> Writing results")
        self.fsi.write()
