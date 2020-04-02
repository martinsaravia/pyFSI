# Open the case file and create the directory structure
#
# Note: the case directory and the case file must have the same name

import importlib, json, pathlib, numpy as np

# Invoke the chosen solver and write run the case
def runCase(caseName):

    # Make the case dictionary
    caseDict = makeCase(caseName)

    # Add the time dict to the case dict
    makeTime(caseDict['control']['time'])

    # Import the solver module
    solverModule = importlib.import_module(
            'solvers.' + caseDict['control']['solver']
            )

    # Call the chosen solver main function (called 'solve')
    solution = getattr(solverModule, 'solve')(caseDict)

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


def makeTime(time):

    # Add the current time and the load factor (f) to the time dict
    time['t'] = np.arange(time['startTime'], time['endTime'], time['deltaT'])
    time['f'] = (time['t'] - time['t'][0]) / time['t'][-1]
