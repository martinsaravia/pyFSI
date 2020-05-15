from properties.materialProperties import constants


def makeDimensionlessNumbers(solid=None, flow=None, fsi=None, numbers=None):
    # Add solid dimensionless numbers
    if numbers is None:  # Mutable object cannot be used as default argument
        numbers = {}     # otherwise a new one is not created (do not use number=None as argument)..
    if flow:
        numbers["Re"] = ReynoldsNumber(flow)
        numbers["Rd"] = ReynoldsNumber(flow, inlet=True)
        numbers["Fr"] = FroudeNumber(flow)
    # Add solid dimensionless numbers
    if solid:

        numbers["Dp"] = displacementNumber(solid)
        numbers["Eg"] = elastoGravityNumber(solid)
    # Add solid dimensionless numbers
    if fsi:
        numbers["Cy"] = CauchyNumber(fsi)
        numbers["Ms"] = massNumber(fsi)
        numbers["Vr"] = reducedVelocityNumber(fsi)

    return numbers


# Dimensionles numbers classes
class dimensionlessNumber(float):
    def __new__(cls, value):
        obj = super().__new__(cls, value)
        return obj

    def info(self):
        print("--> " + self.type + " number is: " + "%10.3E" % self)

    # Base operators
    # def __mul__(self, other):
    #     return self.value * other
    #
    # def __truediv__(self, other):
    #     return self.value / other
    #
    # def __pow__(self, other):
    #     return self.value ** other

# Classes for the different type of dimensionless numbers
# class flowNumber(dimensionlessNumber):
#     def __init__(self):
#         super().__init__()
#
# class solidNumber(dimensionlessNumber):
#     def __init__(self):
#         super().__init__()
#
# class fsiNumber(dimensionlessNumber):
#     def __init__(self):
#         super().__init__()

# Fluid numbers
class ReynoldsNumber(dimensionlessNumber):
    def __new__(cls, flow, inlet=False):
        fluid = flow.fluid()
        if inlet:
            lRef = flow.dRef  # Use the inlet size as reference length
        else:
            lRef = flow.lRef  # Use the flow length as reference length
        value = fluid["rho"] * flow.vRef * lRef / fluid["mu"]
        return super().__new__(cls, value)

    def __init__(self, flow, inlet=False):
        self.type = "Re"

    # Boundary Layer info
    def bLayer(self, geom, dist):
        blthick = 0.0
        blthick = 0.0
        if geom == "plate":
            if self.flowType() == "laminar":
                blthick = 4.91 * dist / self.value**0.5
            elif self.flowType() == "turbulent":
                blthick = 0.37 * dist / self.value**0.2
        return blthick

    # Frow regime info
    def flowRegime(self):
        flowRegime = "laminar"
        if self.value >= 2000:
            flowRegime = "turbulent"
        return flowRegime

class FroudeNumber(dimensionlessNumber):
    def __new__(cls, flow):
        value = flow.vRef / (constants["g"] * flow.lRef) ** 0.5
        return super().__new__(cls, value)

    def __init__(self, flow):
        self.type = "Fr"


# Solid dimensionless numbers
class displacementNumber(dimensionlessNumber):
    def __new__(cls, solid):
        value = solid.uRef / solid.lRef
        return super().__new__(cls, value)

    def __init__(self, solid):
        self.type = "Dp"


class elastoGravityNumber(dimensionlessNumber):
    def __new__(cls, solid):
        mat = solid.material()
        value = mat['rho'] * constants['g'] * solid.lRef / mat['E']
        return super().__new__(cls, value)

    def __init__(self, solid):
        self.type = "Eg"


# FSI numbers
class massNumber(dimensionlessNumber):
    def __new__(cls, fsi):
        value = fsi.flow().fluid()["rho"] / fsi.solid().material()["rho"]
        return super().__new__(cls, value)

    def __init__(self, fsi):
        self.type = "Ms"


class reducedVelocityNumber(dimensionlessNumber):
    def __new__(cls, fsi):
        value = fsi.solid().vRef / fsi.flow().vRef
        return super().__new__(cls, value)

    def __init__(self, fsi):
        self.type = "Vr"


class CauchyNumber(dimensionlessNumber):
    def __new__(cls, fsi):
        value = fsi.flow().fluid()["rho"] * fsi.flow()["vRef"]**2 / fsi.solid().material()["E"]
        return super().__new__(cls, value)

    def __init__(self, fsi):
        self.type = "Cy"

