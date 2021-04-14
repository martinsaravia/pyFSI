# --------------------------------------------------------------------------- #
#    p    #     version: 0.2.0
#    y    #     date: 10/4/2021
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

from copy import deepcopy
from pyFSI.execution.utilities import banner
from pyFSI.execution import io, time
from pyFSI.mesh.fsiMesh1D import fsiMesh1D


class MFSICase:
    def __init__(self, caseName):
        banner("pyFSI v0.3.0 - Martin Saravia 10-3-2021")        # Timing starts

        # Public attributes
        self.DICT = None  # The input dictionary
        self.IODB = None  # The Input-Output Data Base
        self.OREG = None  # The object registry
        self.SOLU = None  # The solver object (yes, weird name)
        self.MFSI = None  # The Magnetic-Fluid-Structure-Interaction object
        self.TIME = None  # The time database

        # Get the path to the caller location
        casePath = pathlib.Path(os.path.abspath(caseName).rsplit('/', 1)[0])  # split cut the name of the file from the Path (added by os)
        print("Running pyFSI case from: ", casePath)

        # Clean the directory tree if present.
        self.clean(casePath)

        # Read the case file and create the fsi object list from the original dictionary (input file)
        self._readInput(casePath, caseName)

        # Create all objects
        self._createObjects()

    def solve(self):
        # Create the case list (contains solver objects, which have a solution)
        print("Calling the solver...")

        self.SOLU.solve()
    
        # End run
        print(f"Problem solved in {self.TIME.elapsed()} seconds ")
        print("************************************************************************")

    # Read the case file
    def _readInput(self, casePath, caseName):
        # casePath = pathlib.Path(__file__).parent.absolute() / ".." / "cases" / caseName
        fileName = caseName + ".json"
        caseFile = casePath / fileName
    
        with open(caseFile, 'r') as f:
            self.DICT = json.load(f)  # Load the input file
    
        # Add paths to the execution dictionary
        paths = {}
        paths['caseName'] = caseName
        paths['casePath'] = casePath
        paths['solidPath'] = casePath / "solid"
        paths['flowPath'] = casePath / "flow"
        paths['fsiPath'] = casePath / "fsi"
        # Copy the solid.0 or fluid.0 fsi.0 to solid or fluid or fsi
        if os.path.isdir(paths['casePath'] / "solid.0"):
            shutil.copytree(paths['casePath'] / "solid.0", paths['solidPath'])
        else:
            paths['solidPath'].mkdir(parents=True, exist_ok=True)
        if os.path.isdir(paths['casePath'] / "flow.0"):
            shutil.copytree(paths['casePath'] / "flow.0", paths['flowPath'])
        else:
            paths['flowPath'].mkdir(parents=True, exist_ok=True)
        if os.path.isdir(paths['casePath'] / "fsi.0"):
            shutil.copytree(paths['casePath'] / "fsi.0", paths['fsiPath'])
        else:
            paths['fsiPath'].mkdir(parents=True, exist_ok=True)

        self.DICT['execution']['paths'] = paths

    def clean(self, casePath, solid=True, flow=True, fsi=True):
        # Get the path to the caller location
        print("Cleaning the case at:", casePath)
        if solid:
            shutil.rmtree(casePath / "solid", ignore_errors=True)
        if flow:
            shutil.rmtree(casePath / "flow", ignore_errors=True)
        if fsi:
            shutil.rmtree(casePath / "fsi", ignore_errors=True)

    def _createObjects(self):
        print("Preprocessing...")

        # Create the object registry
        self.OREG = io.ObjectRegistry()

        self.TIME = time.time(self.DICT["execution"])
        self.OREG.append(self.TIME)

        # Create the mesh object
        mesh = getattr(fsiMesh1D, self.DICT["mesh"]["type"])(self.DICT["mesh"])
        self.OREG.append(mesh)
    
        # Create the solids
        if "solid" in self.DICT:
            print("---> Creating the solid...")
            # Create the solid model
            solidModule = importlib.import_module('pyFSI.models.solidModels.' +
                                                  self.DICT['solid']['formulation'] +
                                                  'Model')
            solid = getattr(solidModule,
                            self.DICT['solid']['formulation'])(self.DICT['execution'],
                                                               self.DICT['solid'],
                                                               mesh,
                                                               self.TIME)
            self.OREG.append(solid)
    
        # Create the boundaries
        if "boundary" in self.DICT:
            boundary = {}
            print("---> Creating the boundaries...")
            # Create the dict of boundary objects
            for b in self.DICT["boundary"]:
                # Import the module (file .py)
                module = importlib.import_module('pyFSI.mesh.region.' + b['type'])
                if "method" in b:
                    # get the class name inside the module
                    obj = getattr(getattr(module, b['type']), b['method'])(mesh, b)
                    # Create the object and store in in the boundary dictionary
                else:
                    # get the class name inside the module
                    obj = getattr(module, b['type'])(solid, b)
    
                self.OREG.append(obj)
    
                boundary[b['name']] = obj
    
        # Create the flows
        if "flow" in self.DICT:
            print("---> Creating the flow...")
            # Create the flow object
            flowModule = importlib.import_module('pyFSI.models.flowModels.' +
                                                 self.DICT['flow']['formulation'] +
                                                 'Model')
            flow = getattr(flowModule,
                           self.DICT['flow']['formulation'])(self.DICT['execution'],
                                                             self.DICT['flow'],
                                                             mesh,
                                                             boundary,
                                                             self.TIME)
            self.OREG.append(flow)
    
        # Create a list of fsi objects
        if "fsi" in self.DICT:
            print("---> Building the FSI object...")
            fsiModule = importlib.import_module('pyFSI.models.fsiModels.' +
                                                self.DICT['fsi']['formulation'] +
                                                'Model')
            mfsi = getattr(fsiModule,
                           self.DICT['fsi']['formulation'])(self.DICT['execution'],
                                                            self.DICT['fsi'],
                                                            solid,
                                                            flow,
                                                            self.TIME)

            self.MFSI = mfsi
            self.OREG.append(mfsi)

        # Create the output data base
        self.IODB = io.IODataBase(self.DICT, self.OREG)
    
        # Create the solver object
        print("---> Building the solver object...")
        solverType = self.DICT['execution']['solver']['type']
        solverModule = importlib.import_module('pyFSI.solvers.' + solverType + 'Solver')
        solver = getattr(solverModule, solverType)(self.MFSI, self.IODB)
        self.SOLU = solver
    
        print("Preprocessing finished successfully...", "\n")
