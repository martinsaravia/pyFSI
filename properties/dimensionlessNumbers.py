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
class dimensionlessNumber:
    def __init__(self):
        pass

    def info(self):
        print("--> " + self.type + " number is: " + "%10.3E" % self.value)

# Fluid numbers
class ReynoldsNumber(dimensionlessNumber):
    def __init__(self, flow, inlet=False):
        super().__init__()
        self.type = "Re"
        fluid = flow.fluid()
        if inlet:
            lRef = flow.dRef  # Use the inlet size as reference length
        else:
            lRef = flow.lRef  # Use the flow length as reference length
        self.value = fluid["rho"] * flow.vRef * lRef / fluid["mu"]

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
        self.value = (fsi.flow().fluid()["rho"] * fsi.flow()["vRef"] ** 2
                      / fsi.solid().material()["E"])

