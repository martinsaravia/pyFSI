import os, importlib, subprocess, numpy as np
from pathlib import Path

from models.solidModels.solidBase import *
from vectors.eigen import eigenValue as eval, eigenVector as evec, eigenSystemVector as esys

import precice as prc

class preciceSolid(object):

    def __repr__(self):
        return 'preciceSolid'

    def __init__(self, execution, control, mesh,  name='NN'):
        self._name = name
        self._fileName = control['filename']
        self.run()

    def run(self):
        # cmd0 = "cd solid;/opt/CalculiX/calculix-adapter/bin/ccx_preCICE -i " + self._fileName + " -precice-participant STRUCTURE;bash"
        # cmd = 'gnome-terminal -- bash -c '  + "'" +  cmd0 + "'"

        cmd = ("gnome-terminal -- bash -c " +
               "'cd solid;" +
               "/opt/CalculiX/calculix-adapter/bin/ccx_preCICE -i " + self._fileName + " -precice-participant STRUCTURE;"
               "bash'")
        # Open a shell and run precice with calculix
        subprocess.call(cmd, shell=True)

