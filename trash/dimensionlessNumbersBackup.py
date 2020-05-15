from properties.materialProperties import constants


def makeDimensionlessNumbers(solid=None, flow=None, fsi=None, numbers={}):
    # Add solid dimensionless numbers
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
        print("--> " + self.type + " number is: " + "%10.3E" % self.value)


# Fluid numbers
class ReynoldsNumber(dimensionlessNumber):
    def __init__(self, flow, inlet=False):
        self.type = "Re"
        fluid = flow.fluid()
        super().__init__()
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
        self.type = "Fr"
        super().__init__()
        self.value = flow.vRef / (constants["g"] * flow.lRef)**0.5


# Solid dimensionless numbers
class displacementNumber(dimensionlessNumber):
    def __new__(cls, solid):
        value = solid.uRef / solid.lRef
        return super().__new__(cls, value)

    def __init__(self, solid):
        self.type = "Dp"


class elastoGravityNumber(dimensionlessNumber):
    def __init__(self, solid):
        self.type = "Eg"
        mat = solid.material()
        super().__init__()
        self.value = mat['rho'] * constants['g'] * solid.lRef / mat['E']

# FSI numbers
class massNumber(dimensionlessNumber):
    def __init__(self, fsi):
        solid = fsi.solid()
        flow = fsi.flow()
        self.type = "Ms"
        super().__init__()
        self.value = flow.fluid()["rho"] / solid.material()["rho"]

class reducedVelocityNumber(dimensionlessNumber):
    def __init__(self, fsi):
        solid = fsi.solid()
        flow = fsi.flow()
        self.type = "Vr"
        super().__init__()
        self.value = solid.vRef / flow.vRef

class CauchyNumber(dimensionlessNumber):
    def __init__(self, fsi):
        solid = fsi.solid()
        flow = fsi.flow()
        self.type = "Cy"
        super().__init__()
        self.value = self.fluid["rho"] * self.fluid["vRef"]**2 / self.solid["E"]
