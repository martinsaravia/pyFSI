import numpy as np, matplotlib.pyplot as plt

#import pycalculix as pyc
from cfd.dimensionlessNumbers.dimensionlessNumber import *
from mesh.fsiMesh1D import fsiMesh1D
from mesh.region.boundary1DRigid import *
from mesh.region.boundary1DBeam import *

from models.solidModels.bernoulliEulerBeamModel.bernoulliEulerBeam import *
from models.solidModels.calculixSolidModel.calculixSolidModel import *
from models.fluidModels.leakageFlow2DModel.leakageFlow2D import *
from models.fsiModels.LFB1dTosiModelMartin2 import LFB1dTosi

# Main function of the solver
def solve(caseDict):

    # Create the mesh
    mesh = getattr(fsiMesh1D, caseDict["mesh"]["type"])(caseDict["mesh"])

    # Get the control dict
    control = caseDict['control']

    # Get the time control from the control dictionary
    time = control['time']

    # Create the solid model
    solid = beamModel(
                control,
                mesh,
                caseDict['solid']
                )

    # Create the boundary dict
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

    # Create the list of fsi objects
    fsi = []
    if caseDict['flow']['bc']['type'] == 'variableInletFlowRate':
        inlet = caseDict['flow']['bc']['inlet']
        # Flow rate as a function of pseudo-time
        for t in time['t']:
            #print("--> Running for time: " + str(t))
            # Modify the inlet velocity before creating the object
            time['ti'] = t
            inlet['Qt'] = inlet['Qi'] + time['ti'] * (inlet['Qf'] - inlet['Qi'])
            # Create the flow object
            flow = leakageFlow2D(
                    control,
                    mesh,
                    boundary,
                    caseDict['flow']
                    )
            # Create the FSI object
            fsi.append(
                LFB1dTosi(
                    control,
                    solid,
                    flow
                    )
                )

    return fsi
