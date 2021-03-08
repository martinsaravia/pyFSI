import importlib
import numpy as np
import os
import subprocess
from pathlib import Path

from pyFSI.models.solidModels.solidBase import *
from pyFSI.vectors.eigen import eigenValue as eval, eigenVector as evec, eigenSystemVector as esys


# import pycalculix as pyc

class calculixBeam(object):

    def __repr__(self):
        return 'calculixBeam '

    def __init__(self, control, mesh, dict, name='NN'):
        self._name = name
        self._control = control
        self._mesh = mesh
        self._dict = dict

        # Init file parameters
        self.fileName = None
        self.fileExt = None
        self.filePath = None
        self._result = None

        # Write the Calculix input file
        self.writeInput(control, mesh, dict)

        # Run Calculix
        self.run()

    def writeInput(self, control, mesh, dict):
        self.fileName = 'solid'
        self.fileExt = self.fileName + '.inp'
        self.filePath = control['paths']['solidPath']
        with open(self.filePath / self.fileExt, 'w') as f:

            # Write node info
            f.write('*NODE, NSET=ALL\n')
            for i, x in enumerate(mesh.x()):
                f.write('\t{0:d}, {1:e}, {2:e}, {3:e}'.format(i+1, x, 0.0, 0.0) + '\n')

            # Write element info
            f.write('*ELEMENT, TYPE=B21, ELSET=' + dict['name'] + '\n')
            for i in range(1, mesh.dim(), 1):
                f.write('\t{0:d}, {1:d}, {2:d}'.format(i, i, i+1) + '\n')

            # Write material info
            dMat = dict['material']
            elasticPars = str(dMat['E']) + ', ' + str(dMat['nu']) + '\n'
            densityPars = str(dMat['rho']) + '\n'
            f.write('*MATERIAL, ' + 'NAME=' + dMat['name'] + '\n')
            f.write('*ELASTIC, TYPE=ISOTROPIC\n')
            f.write('\t' + elasticPars)
            f.write('*DENSITY\n')
            f.write('\t' + densityPars)

            # Write section info
            dSec = dict['section']
            sectionType = dSec['type']
            sectionPars = str(dSec['b']) + ', ' + str(dSec['h'])
            f.write('*BEAM SECTION, ' + 'ELSET=' + dict['name'] + ', ' + 'MATERIAL=' + dMat['name'] + ', ' + 'SECTION=' + sectionType + '\n')
            f.write('\t' + sectionPars + '\n')
            f.write('\t' + '0.0, 1.0, 0.0\n')

            # Write the BC
            f.write('*BOUNDARY\n')
            if dict['bc']['type'] == 'clampedFree':
                f.write('\t1, 1, 6, 0.0\n')
                # f.write('\tALL,3, 3, 0.0\n')  # No es necesaria si uso B21 (2D)

            # Write step InfoArray
            dSol = dict['solution']
            f.write('*STEP\n')
            if dict['solution']['type'] == 'modal':
                f.write('*FREQUENCY, EIGENSOLVER=SUBSPACE\n')
                f.write('\t' + str(dSol['modes']) + ',,,,\n')

            f.write('*NODE FILE, OUTPUT=2D\n')  # Output 2D transforms the beam into a solid.
            f.write('\tU\n')
            f.write('*EL FILE, OUTPUT=2D\n')
            f.write('\tE, S\n')

            # end
            f.write('*END STEP')

    def run(self):
        # runstr = 'ccx ' + str(self.filePath/self.fileName)
        # subprocess.run(runstr, shell=True)
        # print("--> Calculix solver has run correctly...")

        model = pyc.FeaModel(str(self.filePath))  # make the pyCalculix Model

        problem = pyc.Problem(model, 'modal', fname=str(self.filePath/self.fileName)) # Make the problem object

        problem.solveModal()  # Invoke a modified version of pyCalculix

        self._result = problem.rfile.resultsDict()  # Store the results

    def eigenSystem(self):
        # Get the values and vectors from the tricky results structures
        values = list(self._result.keys()) # Eigenvalues as a list of floats
        # print(values)
        # ux = np.zeros(len(mode['node']))
        # for i in mode['node']:
        #     ux[i-1] = mode['node'][i]['ux']
        #     uy[i-1] = mode['node'][i]['uy']
        #     uz[i-1] = mode['node'][i]['uz']
        vectors = []
        theList = list(self._result.values())
        for mode in theList:
            uy = np.zeros(len(mode['node']))
            for i in mode['node']:
                uy[i-1] = mode['node'][i]['uz']  # assume the transverse direction is y
            vectors.append(uy)

        # Eigenvalues and eigenvectors
        evalues = []
        evectors = []
        for i in range(0, self._dict['solution']['modes']):
            # Fill the eigensystem
            val = values[i]
            evalues.append(eval.eigenValue(val))
            # Eigenvectors
            evectors.append(
                evec.eigenVector(
                    vectors[i],
                    self._mesh,
                    info="Beam eigenvector " + str(i+1),
                    normalize=True)
                    )

        eigenSystem = esys.eigenSystemVector(evalues, evectors, self._name)

        return eigenSystem
