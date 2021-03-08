import importlib
import json
import numpy as np
import pathlib
import shutil
import sys
import time as tt
from copy import deepcopy
from functools import partial
from pyFSI.mesh.fsiMesh1D import fsiMesh1D


# Create a list of case dictionaries
def createFSIObjects(caseDict):
    print("Preprocessing...")
    time = caseDict["execution"]["time"]
    caseDictList = []  # List of dictionaries, one for each parameter configuration
    # Process parametric options
    if time["stepping"] == "parametric":
        print("Assembling the parameters vector...")
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
            print("---> Creating the solid...")
            # Create the solid model
            solidModule = importlib.import_module('pyFSI.models.solidModels.' +
                                                  caseDict['solid']['formulation'] +
                                                  'Model')
            solid = getattr(solidModule,
                            caseDict['solid']['formulation'])(caseDict['execution'],
                                                              caseDict['solid'],
                                                              mesh)
        # Create the boundaries
        if "boundary" in caseDict:
            boundary = {}
            print("---> Creating the boundaries...")
            # Create the dict of boundary objects
            for b in caseDict["boundary"]:
                # Import the module (file .py)
                module = importlib.import_module('pyFSI.mesh.region.' + b['type'])
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
            print("---> Creating the flow...")
            # Create the flow object
            flowModule = importlib.import_module('pyFSI.models.flowModels.' +
                                                 caseDict['flow']['formulation'] +
                                                 'Model')
            flow = getattr(flowModule,
                           caseDict['flow']['formulation'])(caseDict['execution'],
                                                            caseDict['flow'],
                                                            mesh,
                                                            boundary)

        # Create a list of fsi objects
        if "fsi" in caseDict:
            print("---> Building the FSI object...")
            fsiModule = importlib.import_module('pyFSI.models.fsiModels.' +
                                                caseDict['fsi']['formulation'] +
                                                'Model')
            fsiObject = getattr(fsiModule,
                                caseDict['fsi']['formulation'])(caseDict['execution'],
                                                                caseDict['fsi'],
                                                                solid,
                                                                flow)

            fsi.append(fsiObject)
            #fsiObject.finish() # Close the files to avoid maximum file error (not good since finish may have other code)

    print("Preprocessing finished successfully...", "\n")
    return fsi