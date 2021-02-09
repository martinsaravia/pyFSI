from models.properties.materialProperties import constants
from models.properties.dimensionlessNumbers import ReynoldsNumber
from tabulate import tabulate
import numpy as np

class boundaryLayer:
    def __init__(self, flow, xPosition=1.0):

        # Initialize attributes
        self.Rex = ReynoldsNumber(flow, L=xPosition)  # x Position: distance from the start of the BL

        # Thickness of the BL
        self.delta = self.thickness()

        # Skin friction coefficient
        self.Cf = self.skinFriction()  # Skin friction

        self.Sw = self.shearStress(flow)

        # Compute the wall distance for yPlus = 1
        self.wD = self.wallDistance(flow)

    def info(self):
        data = [["Reynolds x", self.Rex.value],
                ["Thickness", self.delta],
                ["Friction", self.Cf],
                ["Shear Stress", self.Sw],
                ["Wall y", self.wD]]

        print(tabulate(data,
                       headers=['BL Parameter', 'Value'],
                       tablefmt="fancy_grid",
                       numalign="right"))

    # Wall distance calculation
    def wallDistance(self, flow, yPlus=1):
        # Friction velocity
        fV = (self.Sw / flow.fluid()['rho']) ** 0.5

        # Wall distance
        return yPlus * flow.fluid()['mu'] / (flow.fluid()['rho'] * fV)

    # Calculation of layer thickness
    def thickness(self, geom="plate", formulation="Blasius"):
        delta = None
        if geom == "plate":
            if self.Rex.value <= 1.0E5:  # Laminar boundary layer
                delta = 4.91 * self.Rex.L / self.Rex.value ** 0.5
            else:  # Turbulent boundary layer
                delta = 0.37 * self.Rex.L / self.Rex.value ** 0.2

        return delta

    # Skin friction
    def skinFriction(self):

        # Check if the model is valid
        if self.Rex.value >= 1E9:
            print("--> WARNING! The BL Reynolds number has exceeded "
                  "the max value allowed by the skin friction model")

        Cf = (2 * np.log10(self.Rex.value) - 0.65) ** -2.3

        return Cf

    # Wall shear stress
    def shearStress(self, flow):

        return self.Cf * 0.5 * flow.fluid()['rho'] * self.Rex.V ** 2
