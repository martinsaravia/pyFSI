# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: FSI Region class
#    I    #     return: region
# --------------------------------------------------------------------------- #
# Notes:
#   This class generates a region from two 1D boundary objects
#   Only one boundary can be flexible
#
# --------------------------------------------------------------------------- #

import numpy as np
import scipy.integrate as si



class fsiRegion1D(object):

    def __repr__(self):
        return 'fsiRegion1D'

    def __init__(self, control, mesh, boundary):

        # ----- Public attributes ----- #
        self.name = control['name']
        print(control)
        self.type = control['type']

        # General container for derivatives, integrals, sizes
        self.data = {'s':       np.empty(mesh.size),
                     'six':     np.empty(mesh.size),  # Size indefinite integral
                     'siL':     0,  # Size definite integral between 0-L
                     'ds':      np.empty(mesh.size),  # Vel of the size
                     'dds':     np.empty(mesh.size),  # Vel of the size
                     'dsi':     np.empty(mesh.size),  # Vel of the size integral
                     'ddsi':    np.empty(mesh.size)  # Accel of the size integral
                     }

        # ----- Private attributes ----- #
        self._mesh = mesh   # Reference to the mesh
        self._bBot = boundary[control['botBoundary']]  # Top boundary
        self._bTop = boundary[control['topBoundary']]   # Bottom boundary

        # ----- Procedures ----- #
        # Find the associated flexible boundary
        if self._bBot.isFlexible():
            flexBoundary = self._bBot
        elif self._bTop.isFlexible():
            flexBoundary = self._bTop
        else:
            flexBoundary = None
            print("---> Warning: No flexible boundary found for region" + self.name)

        self._eigen = flexBoundary.eigen()

        self.check()

    # Update the geometric data
    def update(self):
        self._bTop.update()
        self._bBot.update()
        self.data['s'] = self._bTop.y - self._bBot.y
        self.data['six'] = self._bTop.ix - self._bBot.ix
        self.data['siL'] = self._bTop.iL - self._bBot.iL
        self.data['ds'] = self._bTop.dy - self._bBot.dy
        self.data['dds'] = self._bTop.ddy - self._bBot.ddy
        self.data['dsi'] = self._bTop.dyi - self._bBot.dyi
        self.data['ddsi'] = self._bTop.ddyi - self._bBot.ddyi



    def check(self):
        # Check if the mesh density is ok
        if (self._mesh.x[1] - self._mesh.x[0]) > self.data['s'][0]:
            print("--> WARNING: mesh dx is smaller the channel inlet "
                  "for region " + self.name + ". Results may be wrong. ")

    # ----- Getters ----- #
    # Return reference to the mesh
    def mesh(self):
        return self._mesh

    # Return reference to the flex boundary (only 1 flex boundary)
    def eigen(self):
        return self._eigen

    # Reference to the boundary objects
    def top(self):
        return self._bTop

    def bot(self):
        return self._bBot


