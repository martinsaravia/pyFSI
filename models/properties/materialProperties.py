
fluids = {
    "air": {
        "rho": 1.204,
        "mu":  1.825E-5,
        "nu":  1.516E-5
    },

    "water": {
        "rho": 0.99802E+3,
        "mu":  1.002E-3,
        "nu":  1.004E-6
    }
}


solids = {

    "steel": {
        "rho": 7.85E+3,
        "nu":  0.3,
        "E":   2.1E11
    },

    "aluminum": {
        "rho": 2.7E+3,
        "nu":  0.33,
        "E":   7.2E+10
    },

    "kaneko": {
        "rho": 8.78E+3,
        "nu":  0.3,
        "E":   1.1E+11
    }

}


constants = {
    "g": 9.8
}


def convert(solidDict=None, fluidDict=None, constDict=None, unit="m"):
    # Unit conversion
    if unit == "mm":
        print("--> Converting constant units to mm.")
        if constDict:
            constDict["g"] *= 9800.0
        if fluidDict:
            print("--> Converting fluid units to mm.")
            fluidDict["rho"] *= 1.0E-9
            fluidDict["mu"] *= 1.0E-3
        if solidDict:
            print("--> Converting solid units to mm.")
            solidDict["rho"] *= 1.0E-9
            solidDict["E"] *= 1.0E-3
