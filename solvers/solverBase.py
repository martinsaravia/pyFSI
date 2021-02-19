# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 07/12/2020
#    F    #     author: Martin Saravia
#    S    #     description: Base class for all solvers
#    I    #
# --------------------------------------------------------------------------- #

import importlib
from abc import ABCMeta, abstractmethod

from mesh.fsiMesh1D import fsiMesh1D


class solverBase:
    def __init__(self, fsi):
        # Public Attributes
        self.fsi = fsi
        self.time = None
        self.control = fsi.execution()

        # Private attributes
        # self._ObjectRegistry = []  # The object registry
        # self._caseDict = caseDict
        #
        # self.createObjects(caseDict)

    # Create the FSI objects (meshes, solids, flows, boundaries)
    def createObjects(self, caseDict):
        # Create the mesh
        print("Creating the meshes...")
        # Only one mesh type is available.
        mesh = getattr(fsiMesh1D, caseDict["mesh"]["type"])(caseDict["mesh"])
        self._ObjectRegistry.append(mesh)  # Add the mesh to the object registry

        # Create the solids
        if "solid" in caseDict:
            print("Creating the solids...")
            # Create the solid model
            solidModule = importlib.import_module('models.solidModels.' +
                                                  caseDict['solid']['formulation'] +
                                                  'Model')
            solid = getattr(solidModule,
                            caseDict['solid']['formulation'])(caseDict['execution'],
                                                              caseDict['solid'],
                                                              mesh)
            self._ObjectRegistry.append(solid)  # Add the solid to the object registry

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

                self._ObjectRegistry.append(obj)  # Add the boundary to the object registry

        # Create the flows
        if "flow" in caseDict:
            print("Creating the flows...")
            # Create the flow object
            flowModule = importlib.import_module('models.flowModels.' +
                                                 caseDict['flow']['formulation'] +
                                                 'Model')
            flow = getattr(flowModule,
                                caseDict['flow']['formulation'])(caseDict['execution'],
                                                                 caseDict['flow'],
                                                                 mesh,
                                                                 boundary)
            self._ObjectRegistry.append(flow)  # Add the flow to the object registry

        # Create a list of fsi objects
        if "fsi" in caseDict:
            print("Building the FSI models...")
            fsiModule = importlib.import_module('models.fsiModels.' +
                                                caseDict['fsi']['formulation'] +
                                                'Model')
            self.fsi.append(getattr(fsiModule,
                                    caseDict['fsi']['formulation'])(caseDict['execution'],
                                                                    caseDict['fsi'],
                                                                    solid,
                                                                    flow))

    # Execute final tasks for the fsi objects
    def finish(self):
        self.fsi.finish()

    # Access functions
    def caseDict(self):
        return self._caseDict

    # Abstract methods
    @abstractmethod
    def solve(self):
        pass

    @abstractmethod
    def advance(self, tspan):
        pass

    # Write the output of every object
    @abstractmethod
    def write(self):
        pass

