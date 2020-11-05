# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: BE Beam class
#    I    #     return: solidModel object
# --------------------------------------------------------------------------- #
# Notes:
#   Staggered solver for fsi
#   The solid integrated with RK45 with time subdivition is dissipative but if
#   we integrate the whole time with RK45 it is conservative. It appears that
#   there is problem with the update of the initial condition but Radau has not
#   this problem, so I may well be a problem of the algorithm.
# --------------------------------------------------------------------------- #
# Change the system matrix.
# Corrected the position of G, Small acceleration and changed the sign of Tb

# import pycalculix as pyc
import importlib
from mesh.fsiMesh1D import fsiMesh1D
import scipy.integrate as si
import matplotlib.pyplot as plt
import numpy as np


# from models.solidModels import bernoulliEulerBeamModel


# Main function of the solver
def solve(caseDict):
    # Create the mesh
    mesh = getattr(fsiMesh1D, caseDict["mesh"]["type"])(caseDict["mesh"])

    # Create the solid model
    solidModule = importlib.import_module('models.solidModels.' +
                                          caseDict['solid']['formulation'] +
                                          'Model')

    solid = getattr(solidModule,
                    caseDict['solid']['formulation'])(
        caseDict['execution'],
        caseDict['solid'],
        mesh)

    # Create the dict of boundary objects
    boundary = {}
    for b in caseDict["boundary"]:
        # Import the module (file .py)
        module = importlib.import_module('mesh.region.' + b['type'])
        if "method" in b:
            # get the class name inside the module
            cls = getattr(getattr(module, b['type']), b['method'])(mesh, b)
            # Create the object and store in in the boundary dictionary
            boundary[b['name']] = cls
        else:
            # get the class name inside the module
            cls = getattr(module, b['type'])(solid, b)
            boundary[b['name']] = cls

    # Create the flow object
    flowModule = importlib.import_module('models.flowModels.' +
                                         caseDict['flow']['formulation'] +
                                         'Model')
    flow = getattr(flowModule,
                   caseDict['flow']['formulation'])(caseDict['execution'],
                                                    caseDict['flow'],
                                                    mesh,
                                                    boundary)

    # Create a list of fsi objects
    fsi = []
    fsiModule = importlib.import_module('models.fsiModels.' +
                                        caseDict['fsi']['formulation'] +
                                        'Model')
    fsi.append(getattr(fsiModule,
                       caseDict['fsi']['formulation'])(caseDict['execution'],
                                                       caseDict['fsi'],
                                                       solid,


                                                       flow))
    # Solver  loop
    time = caseDict['execution']['time']
    solver = caseDict['execution']['solver']

    tspan = (caseDict['execution']['time']['startTime'],
             caseDict['execution']['time']['endTime'])

    initialConditions = fsi[0].state
    for ti in range(1, len(time['t'])):
        for i in range(1):
            tspan = [time['t'][ti-1], time['t'][ti]]
            dt = time['t'][ti] - time['t'][ti-1]
            print('---> Solving for time: ', time['t'][ti])

            fSol = si.solve_ivp(fsi[0].flow().rhs,
                                tspan,
                                fsi[0].fState,
                                method='Radau')

            fsi[0].update('fluid', time['t'][ti], fSol.y[:, -1])

            # Solve the solid
            sSol = si.solve_ivp(fsi[0].rhsSolid,
                                tspan,
                                fsi[0].sState,
                                method='Radau')
            fsi[0].update('solid', sSol.t[-1], sSol.y[:, -1])

            # Update the region geometry
            fsi[0].update('regions', None,  None)

            # # Solve the Fluid
            fsi[0].update('bc', time['t'][ti], None)

            fSol = si.solve_ivp(fsi[0].flow().rhs,
                                tspan,
                                fsi[0].fState,
                                method='Radau')

            fsi[0].update('fluid', time['t'][ti], fSol.y[:, -1])

        # fsi[0].residual()

        # plt.figure()
        # plt.plot(fsi[0].flow().px[0])
        # plt.plot(fsi[0].flow().px[1])
        # # plt.figure()
        # print('fSol', fSol.y)
        #
        # plt.show()

        a= 1
        # dp = fsi[0].flow().px[0] - fsi[0].flow().px[1]
        # plt.figure()
        # plt.plot(dp)
        # v = fsi[0].solid().eigen.vectors
        # x = fsi[0].solid().mesh().x
        # # Modal pressure loads
        # lm1 = dp * v[0]
        # # lm2 = dp * v[1]
        # plt.figure()
        # plt.plot(lm1)
        # # plt.plot(lm2)
        #
        # print('i1: ', si.simps(lm1, x))
        # print('i2: ', si.simps(lm2, x))

    fsi[0].finish()

    return fsi
