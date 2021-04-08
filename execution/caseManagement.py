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
import importlib
import json
import pathlib
import shutil
import sys
import os
import time as tt
from copy import deepcopy
from pyFSI.execution.objectCreation import *
from pyFSI.execution.utilities import banner


def runCase(caseName):
    banner("pyFSI v0.2.1 - Martin Saravia 23-2-2021")

    # Timing starts
    start = tt.perf_counter()

    # Get the path to the caller locatiion
    casePath = pathlib.Path(os.path.abspath(caseName).rsplit('/', 1)[0])# Rsplit cut the name of the file from the Path (added by os)
    print("Running pyFSI case from: ", casePath)

   # Read the case file

    caseDict = readCase(casePath, caseName)
    # Create the fsi object list from the original dictionary (input file)
    fsiObjects = createFSIObjects(caseDict)

    # Create the case list (contains solver objects, which have a solution)
    print("Calling the solver...")
    case = []
    for fsi in fsiObjects:
        # Create the solver object and solve the case, the solver has a solution
        solverType = fsi.execution()['solver']['type']
        solverModule = importlib.import_module('pyFSI.solvers.' + solverType + 'Solver')
        solver = getattr(solverModule, solverType)(fsi)
        solver.solve()
        case.append(solver)

    # End run
    finish = tt.perf_counter()
    print(f"Problem solved in {round(finish - start, 2)} seconds ")
    print("************************************************************************")

    return case

# Read the case file
def readCase(casePath, caseName):
    # casePath = pathlib.Path(__file__).parent.absolute() / ".." / "cases" / caseName
    fileName = caseName + ".json"
    caseFile = casePath / fileName

    with open(caseFile, 'r') as f:
        caseDict = json.load(f)  # Load the input file

    # Add paths to the execution dictionary
    paths = {}
    paths['caseName'] = caseName
    paths['casePath'] = casePath
    paths['solidPath'] = casePath / "solid"
    paths['fluidPath'] = casePath / "flow"
    paths['fsiPath'] = casePath / "fsi"
    paths['solidPath'].mkdir(parents=True, exist_ok=True)
    paths['fluidPath'].mkdir(parents=True, exist_ok=True)
    paths['fsiPath'].mkdir(parents=True, exist_ok=True)
    caseDict['execution']['paths'] = paths

    return caseDict

def cleanCase(caseName, solid=True, flow=True, fsi=True):
    # Get the path to the caller locatiion
    casePath = pathlib.Path(os.path.abspath(caseName).rsplit('/', 1)[0])# Rsplit cut the name of the file from the Path (added by os)
    print("Cleaning the case at:", casePath)
    if solid:
        shutil.rmtree(casePath / "solid", ignore_errors=True)
    if flow:
        shutil.rmtree(casePath / "flow", ignore_errors=True)
    if fsi:
        shutil.rmtree(casePath / "fsi", ignore_errors=True)