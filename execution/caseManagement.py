# --------------------------------------------------------------------------- #
#    p    #     version: 0.2.0
#    y    #     date: 10/2/2021
#    F    #     author: Martin Saravia
#    S    #     description: routines for case management
#    I    #     return:
# --------------------------------------------------------------------------- #
# Notes:
# Available functions are readCase, cleanCase and runCase
# --------------------------------------------------------------------------- #
import importlib, json, pathlib, sys, shutil
import time as tt
from copy import deepcopy
from execution.utilities import banner
from execution.objectCreation import *


def runCase(caseName):

    banner("pyFSI v0.2 - Martin Saravia 21-12-2020")

    # Timing starts
    start = tt.perf_counter()

    # Read the case file
    caseDict = readCase(caseName)

    # Create the fsi object list from the original dictionary (input file)
    fsiObjects = createFSIObjects(caseDict)

    # Create the case list (contains solver objects, which have a solution)
    print( "------------------------------------------------------------------------")
    print("Solving...")
    case = []
    for fsi in fsiObjects:
        # Create the solver object and solve the case, the solver has a solution
        solverType = fsi.execution()['solver']['type']
        solverModule = importlib.import_module('solvers.' + solverType + 'Solver')
        solver = getattr(solverModule, solverType)(fsi)
        solver.solve()
        case.append(solver)

    # End run
    finish = tt.perf_counter()
    print(f"Problem solved in {round(finish - start, 2)} seconds ")
    print("************************************************************************")

    return case

# Read the case file
def readCase(caseName):
    casePath = pathlib.Path(__file__).parent.absolute() / ".." / "cases" / caseName
    caseFile = caseName + ".json"
    with open(casePath / caseFile, 'r') as f:
        caseDict = json.load(f)  # Load the input file

    # Add paths to the execution dictionary
    paths = {}
    paths['caseName'] = caseName
    paths['casePath'] = casePath
    paths['solidPath'] = casePath / "solid/"
    paths['fluidPath'] = casePath / "fluid/"
    paths['fsiPath'] = casePath / "fsi/"
    paths['solidPath'].mkdir(parents=True, exist_ok=True)
    paths['fluidPath'].mkdir(parents=True, exist_ok=True)
    paths['fsiPath'].mkdir(parents=True, exist_ok=True)
    caseDict['execution']['paths'] = paths

    return caseDict

def cleanCase():
    shutil.rmtree("solid", ignore_errors=True)
    shutil.rmtree("flow", ignore_errors=True)
    shutil.rmtree("fsi", ignore_errors=True)