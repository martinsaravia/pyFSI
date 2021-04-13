""" Input-Output Module
This module holds the classes used for writing input and output data.
"""
import numpy as np
from pyFSI.execution.errors import error
import sys

class IOFile:
    """
    Class for output files. It takes and object and a variable and creates an output file that is written when the database
    write method is called.
    """
    def __init__(self, obj, variable, mode='a+', bufferSize=1):
        filename = variable + ".out"
        self.location = obj.path / filename
        self.obj = obj  # Object from which we read the data
        try:
            self.var = obj.varMap[variable]
        except KeyError:
            error("ERROR: Variable: " + variable + " is not available for output for object named: " + obj.name + "\n" +
                  "       Valid variables are: " + str(list(obj.varMap.keys()))[:])

        self.expr = "self.obj." + self.var  # Expression to evaluate the variable value

        self.file = open(self.location,
                         mode,
                         buffering=bufferSize)

    def write(self):
        value = eval(self.expr)
        if isinstance(value, np.floating) or isinstance(value, float) or isinstance(value, int):
            self.file.write(str(value) + '\n')
        else:
            self.file.write(" ".join(map(str, value)) + '\n')  # If list, split it

    def close(self):
        self.file.close()


class IODataBase:
    """
    Input-Output data base. I holds all the output files (type IOFile).
    """
    def __init__(self, caseDict, registry):
        self.files = []
        # Add current time time to the database
        time = registry.get('time')
        self.files.append(IOFile(time, 'time'))
        # Add the rest of the variables to the database
        for key, val in caseDict['execution']['output'].items():
            obj = registry.get(key)  # Find the object by name
            for i in val:  # Iterate through variable names of the current object
                file = IOFile(obj, i)
                self.files.append(file)

    def write(self):
        for file in self.files:
            file.write()

    def close(self):
        for file in self.files:
            file.close()

class objectRegistry:
    """
    A class for registering all objects. I has a get method for retrieving an object by name.
    """
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


class Logger:
    """
    Class for logging simultaneously to the console and to a file.
    From https://stackoverflow.com/questions/616645/how-to-duplicate-sys-stdout-to-a-log-file
    """
    def __init__(self, name='stdout.log', mode='w'):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self

    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        self.file.flush()

