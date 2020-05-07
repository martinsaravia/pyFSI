# Open the case file and create the directory structure
#
# Note: the case directory and the case file must have the same name

import importlib, json, pathlib
import numpy as np
from copy import deepcopy
# Invoke the solver and write run the case
def runCase(caseName):

    # Make the main case dictionary structure
    caseDict = makeCase(caseName)

    # Add the time dict to the case dict
    makeTime(caseDict['control'])

    # Import the solver module
    solverModule = importlib.import_module(
            'solvers.' + caseDict['control']['solver']
            )

    solution = []

    for p in caseDict['control']['parameters']['p']:
        # Copy the case dict in order to fill it with parametric stuff
        newCaseDict = deepcopy(caseDict)  # I know, bad practice
        pDict = newCaseDict['control']['parameters']
        pDict['pi'] = p  # Current Parameter value
        print("-----------------------------------------------")
        print("-       Running for parameter: " + str(pDict['pi']) + "          -")
        print("-----------------------------------------------")
        # Parametric boundary
        if pDict['type'] == "boundary":
            for name in pDict['names']:
                for idx, boun in enumerate(newCaseDict['boundary']):
                    if boun['name'] == name:
                        for var in pDict['vars']:
                            boun[var] = pDict['pi'] * newCaseDict['boundary'][idx][var]
            # Call the chosen solver main function (called 'solve') for this parameter configuration
            solution.append(getattr(solverModule, 'solve')(newCaseDict))

    return solution


# Make the case directory tree and create the case dictionary
def makeCase(caseName):
    casePath = pathlib.Path(__file__).parent.absolute()/".."/"cases"/caseName
    caseName += ".json"
    with open(casePath/caseName, 'r') as f:
        caseDict = json.load(f) # Load the input file

    # Add paths to the control dictionary
    paths = {}
    paths['casePath'] = casePath
    paths['solidPath'] = casePath/"solid/"
    paths['fluidPath'] = casePath/"fluid/"
    paths['solidPath'].mkdir(parents=True, exist_ok=True)
    paths['fluidPath'].mkdir(parents=True, exist_ok=True)
    caseDict['control']['paths'] = paths

    return caseDict


def makeTime(control):

    # Add the current time and the load factor (f) to the time dict
    time = control['time']
    time['t'] = np.arange(time['startTime'], time['endTime']+time['deltaT'], time['deltaT'])  # Real Time
    time['steps'] = len(time['t'])
    time['ti'] = time['t'][0]  # Initial time
    
    # Support for parametric studies
    para = control['parameters']
    if 'parameters' in control:
        para['p'] = np.linspace(para['iPF'], para['fPF'], para['steps'])  # Proportional factor vector
    else:
        para['p'] = [1]
    para['pi'] = para['p'][0]  # Initial Parameter
    
