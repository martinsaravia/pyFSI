from solvers.solverBase import solverBase
import scipy.integrate as si
import numpy as np

# Solver for transient FSI simulations
class transient(solverBase):
    def __init__(self, fsi):
        super().__init__(fsi)

    def solve(self):
        time = self.control['time']
        initialConditions = self.fsi.state
        tspan = np.array([time['startTime'], time['startTime'] + time['deltaT']])
        dt = time['deltaT']
        while tspan[1] <= time['endTime']:
            print("-> Solving for time: ", tspan[1])
            self.advance(tspan)
            tspan += dt
        self.finish()

    def advance(self, tspan):
        fsi = self.fsi

        # 1) Solve the flow
        print("---> Solving the flow")
        fSol = si.solve_ivp(fsi.flow().rhs,
                            tspan,
                            fsi.fState,
                            method='Radau')
        fsi.update('flow', tspan[1], fSol.y[:, -1])

        # 2) Update the boundary condition
        print("---> Updating the flow force over the solid")
        fsi.update('forces', tspan[1], None)

        # 3) Solve the solid
        print("---> Solving the solid")
        sSol = si.solve_ivp(fsi.rhsSolid,
                            tspan,
                            fsi.sState,
                            method='Radau')
        fsi.update('solid', tspan[1], sSol.y[:, -1])

        # 4) Update the region geometry
        print("---> Updating the region position, velocities and accelerations")
        fsi.update('regions', None, None)

        self.write(tspan[1])

    def write(self, t):
        print("---> Writing results for time ", t)
        self.fsi.write()
