# --------------------------------------------------------------------------- #
#    p    #     version: 0.2
#    y    #     date: 21/122020
#    F    #     author: Martin Saravia
#    S    #     description: run script
#    I    #     return: solver object
# --------------------------------------------------------------------------- #
# Notes:
#   Staggered solver for fsi using preCICE and Calculix
#   a) The solid is solved using Calculix
# --------------------------------------------------------------------------- #
import importlib, json, pathlib, sys, shutil
from functools import partial
import time as tt
from mesh.fsiMesh1D import fsiMesh1D
import numpy as np
from copy import deepcopy


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

# Nice banner to print at the start (or not to nice)
def banner(s, width=69):
    stars = '*' * width
    pad = (width + len(s)) // 2
    print(f'{stars}\n{s:>{pad}}\n{stars}')


# Create a list of case dictionaries
def createFSIObjects(caseDict):
    time = caseDict["execution"]["time"]
    caseDictList = []  # List of dictionaries, one for each parameter configuration
    # Process parametric options
    if time["stepping"] == "parametric":
        print("This is a parametric analysis...")
        parValues = np.linspace(time["startParameter"],
                                time["endParameter"],
                                time["steps"])
        np.savetxt("fsi/parameter.out", parValues) # Save the parameter to a file
        # Choose the parameter
        for par in parValues:
            newCaseDict = deepcopy(caseDict)  # I know, bad practice
            newCaseDict["fsi"]["parameter"] = par # Store the parameter
            if time["parameter"] == "inletFlowRate":
                # Assign the parameter to the variable.
                newCaseDict["flow"]["bc"]["inlet"]["Q0"] = par
                # Change the name of the FSI object
                newCaseDict["fsi"]["name"] += ". Parameter: " + time["parameter"] + " " + str(par)
            caseDictList.append(newCaseDict)

    else:
        caseDictList.append(caseDict)  # If not parametric just add the case dictionary to make a 1 element list

    # Create all FSI objects
    fsi = []
    for caseDict in caseDictList:
        # Only one mesh type is available.
        mesh = getattr(fsiMesh1D, caseDict["mesh"]["type"])(caseDict["mesh"])

        # Create the solids
        if "solid" in caseDict:
            print("Creating the solid...")
            # Create the solid model
            solidModule = importlib.import_module('models.solidModels.' +
                                                  caseDict['solid']['formulation'] +
                                                  'Model')
            solid = getattr(solidModule,
                            caseDict['solid']['formulation'])(caseDict['execution'],
                                                              caseDict['solid'],
                                                              mesh)
        # Create the boundaries
        if "boundary" in caseDict:
            boundary = {}
            print("Creating the boundaries...")
            # Create the dict of boundary objects
            for b in caseDict["boundary"]:
                # Import the module (file .py)
                module = importlib.import_module('mesh.region.' + b['type'])
                if "method" in b:
                    # get the class name inside the module
                    obj = getattr(getattr(module, b['type']), b['method'])(mesh, b)
                    # Create the object and store in in the boundary dictionary
                else:
                    # get the class name inside the module
                    obj = getattr(module, b['type'])(solid, b)

                boundary[b['name']] = obj

        # Create the flows
        if "flow" in caseDict:
            print("Creating the flow...")
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
        if "fsi" in caseDict:
            print("Building the FSI object...")
            fsiModule = importlib.import_module('models.fsiModels.' +
                                                caseDict['fsi']['formulation'] +
                                                'Model')
            fsi.append(getattr(fsiModule,
                                caseDict['fsi']['formulation'])(caseDict['execution'],
                                                                caseDict['fsi'],
                                                                solid,
                                                                flow))

    print("All objects were created succesfully...", "\n")
    return fsi