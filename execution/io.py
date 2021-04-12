import numpy as np


class IOObject:
    def __init__(self, obj, variable, mode='a+', bufferSize=1):
        filename = variable + ".out"
        self.location = obj.path / filename
        self.obj = obj  # Object from which we read the data
        self.var = obj.varMap[variable]

        self.expr = "self.obj." + self.var  # Expression to evaluate the variable value

        self.file = open(self.location,
                         mode,
                         buffering=bufferSize)

    def write(self):
        value = eval(self.expr)
        if isinstance(value, np.floating) or isinstance(value, float):
            self.file.write(str(value) + '\n')
        else:
            self.file.write(" ".join(map(str, value)) + '\n')  # If list, split it

    def close(self):
        self.file.close()

# Input-Output data base
class IODataBase:
    def __init__(self, caseDict, registry):
        self.files = []
        # Add current time time to the database
        time = registry.get('time')
        self.files.append(IOObject(time, 'value'))
        # Add the rest of the variables to the database
        for key, val in caseDict['execution']['output'].items():
            obj = registry.get(key)  # Find the object by name
            for i in val:  # Iterate through variable names of the current object
                file = IOObject(obj, i)
                self.files.append(file)

    def write(self):
        for file in self.files:
            file.write()

    def close(self):
        for file in self.files:
            file.close()


class objectRegistry:
    def __init__(self):
        self.objects = []
        self.names = []

    def append(self, obj):   # Add an object to the registry
        self.objects.append(obj)
        self.names.append(obj.name)

    def get(self, name):  # Get and object from the registry
        for i in self.objects:
            if i.name == name:
                return i
        # Raise error if object is not found
        err_msg = "Object named: " + name + " not found in list" + str(self.names)
        raise RuntimeError(err_msg)

# Variable mapping
varmap = {
    "value":                "value",
    "flowRates":            "Q0",
    "flowSpeeds":           "v0",
    "naturalFrequencies":   "eigen.values",
    "eigenValues":          "ES.evalues()",
    "eigenVectors":         "ES.evectors()",
    
    "displacements":        "y['mid']",
    "velocities":           "dy['mid']",
    "accelerations":        "ddy['mid']",
    "numbers":              "numbers"
}
