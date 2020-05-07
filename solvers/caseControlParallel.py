# Open the case file and create the directory structure
#
# Note: the case directory and the case file must have the same name

import importlib, json, pathlib
import numpy as np
from copy import deepcopy
from multiprocessing import Pool
from itertools import repeat
import time
from tqdm import tqdm
from functools import partial


# Invoke the solver and write run the case
def runCase(caseName):

    banner("pyFSI v0.1 - Martin Saravia 05-05-2020")


    start = time.perf_counter()

    # Make the main case dictionary structure and return the case Dict
    caseDict = makeCase(caseName)

    # Add the time dict to the case dict
    makeTime(caseDict['control'])

    # Import the solver module
    solverModule = importlib.import_module(
            'solvers.' + caseDict['control']['solver']
            )

    # Create a dictionary per parameter value and store in a list
    caseDictList = []  # List of dictionaries, one for each parameter configuration
    for p in caseDict['control']['parameters']['p']:
        # Copy the case dict in order to fill it with parametric stuff
        newCaseDict = deepcopy(caseDict)  # I know, bad practice
        pDict = newCaseDict['control']['parameters']
        pDict['pi'] = p  # Current Parameter value
        # Parametric boundary
        if pDict['type'] == "boundary":
            for name in pDict['names']:
                for idx, boun in enumerate(newCaseDict['boundary']):
                    if boun['name'] == name:
                        for var in pDict['vars']:
                            boun[var] = pDict['pi'] * newCaseDict['boundary'][idx][var]
        caseDictList.append(newCaseDict)


    # Solve serial or parallel
    # global progressBar
    # progressBar = tqdm(total=len(caseDictList))
    solution = []
    # Start the processing and map to the list of input dictionaries
    if "processes" in caseDict['control']:  # parallel run
        nproc = caseDict['control']["processes"]


        pool = Pool(nproc)

        # Asyncronous variant
        # solution = pool.map_async(parallel, caseDictList).get()

        # Asyncronour unordered
        # res = pool.imap_unordered(parallel, caseDictList)
        # for _ in tqdm(res, total=len(caseDictList)):
        #     pass
        # pool.close()
        # pool.join()
        # solution = res

        # Syncronous variant
        with Pool(nproc) as p:
            solution = p.map(parallel, caseDictList)  # SImple map
            p.close()
            p.join()

    else:
        nproc = 1
        for d in caseDictList:
            solution.append(getattr(solverModule, 'solve')(newCaseDict))

    # Print elapsed time and return solution
    finish = time.perf_counter()
    print(f"\n--> Problem solved in {round(finish - start, 2)} seconds with {nproc} processes")
    return solution


def parallel(iCaseDict):
    solverModule = importlib.import_module(
            'solvers.' + iCaseDict['control']['solver']
            )
    sol = getattr(solverModule, 'solve')(iCaseDict)
    #progressBar.update(8)
    print("--> Solving for parameter: " + str(round(iCaseDict['control']['parameters']['pi'],3)),end='\r')
    return sol


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
    


def banner(s, width=69):
    stars = '*' * width
    pad = (width + len(s)) // 2
    print( f'{stars}\n{s:>{pad}}\n{stars}')