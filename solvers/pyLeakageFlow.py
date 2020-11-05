# import pycalculix as pyc
import importlib
from mesh.fsiMesh1D import fsiMesh1D
#from models.solidModels import bernoulliEulerBeamModel


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

    # Create a list of fsi objects
    fsi = []
    if caseDict['flow']['bc']['type'] == 'variableInletFlowRate':
        inlet = caseDict['flow']['bc']['inlet']
        # Flow rate as a function of pseudo-time
        for t in caseDict['execution']['time']['t']:
            # Modify the inlet velocity before creating the object
            inlet['Qt'] = inlet['Qi'] + t * (inlet['Qf'] - inlet['Qi'])
            caseDict['execution']['time']['ti'] = t

            # Create the flow object
            flowModule = importlib.import_module('models.flowModels.' +
                                                 caseDict['flow']['formulation'] +
                                                 'Model')
            flow = getattr(flowModule,
                           caseDict['flow']['formulation'])(
                           caseDict['execution'],
                           caseDict['flow'],
                           mesh,
                           boundary)

            # Create the FSI object
            fsiModule = importlib.import_module('models.fsiModels.' +
                                                caseDict['fsi']['formulation'] +
                                                'Model')
            fsi.append(getattr(fsiModule,
                               caseDict['fsi']['formulation'])(
                               caseDict['execution'],
                               caseDict['fsi'],
                               solid,
                               flow))

    return fsi
