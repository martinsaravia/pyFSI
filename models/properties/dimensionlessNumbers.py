""" This module contains all the classes for calculating dimensionless numbers"""

from pyFSI.models.properties.materialProperties import constants


class dimensionlessNumber:
    """ Base class for dimensionless numbers"""
    def __init__(self, value=None):
        self.value = value

    def __str__(self):
        return str(self.value)

    def info(self):
        print("--> " + self.type + " number is: " + "%10.3E" % self.value)


class ReynoldsNumber(dimensionlessNumber):
    def __init__(self, flow, L=None, V=None, type="Re"):
        super().__init__()
        self.type = type

        # Get the reference velocity
        if V is not None:
            self.V = V
        else:
            self.V = flow.vRef  # Reference velocity of the flow

        # Choose the reference length (if lRef is given,
        # it takes precedence over type value)
        if L is not None:
            self.L = L
        else:
            if self.type == "Re":  # domain length Reynolds
                self.L = flow.lRef
            elif self.type == "Rd":  # inlet dimension Reynolds
                self.L = flow.dRef

        self.value = flow.fluid()["rho"] * self.V * self.L / flow.fluid()["mu"]

    # Frow regime info
    def flowRegime(self):
        flowRegime = "laminar"
        if self.value >= 2000:
            flowRegime = "turbulent"
        return flowRegime


class FroudeNumber(dimensionlessNumber):
    def __init__(self, flow):
        super().__init__()
        self.type = "Fr"
        self.value = flow.vRef / (constants["g"] * flow.lRef) ** 0.5


# Solid dimensionless numbers
class displacementNumber(dimensionlessNumber):
    def __init__(self, solid):
        super().__init__()
        self.type = "Dp"
        self.value = solid.uRef / solid.lRef

class elastoGravityNumber(dimensionlessNumber):
    def __init__(self, solid):
        super().__init__()
        self.type = "Eg"
        mat = solid.material()
        self.value = mat['rho'] * constants['g'] * solid.lRef / mat['E']


# FSI numbers
class massNumber(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Ms"
        self.value = fsi.flow().fluid()["rho"] / fsi.solid().material()["rho"]

class reducedVelocityNumber(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Vr"
        self.value = fsi.solid().vRef / fsi.flow().vRef

class CauchyNumber(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Cy"
        self.value = (fsi.flow().fluid()["rho"] * fsi.flow().vRef**2 /
                      fsi.solid().material()["E"])

