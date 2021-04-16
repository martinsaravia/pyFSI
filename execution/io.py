""" Input-Output Module
This module holds the classes used for writing input and output data.
"""
import sys
from abc import ABC, abstractmethod
import numpy as np
from pyFSI.execution.errors import error

class IOFile:
    """
    Class for output files. It takes and object and a variable and creates an output file that is written when the database
    write method is called.
    """
    def __init__(self, obj, variable, mode='a+', bufferSize=1):
        filename = variable + ".out"
        self.location = obj.path / filename
        self.obj = obj  # Object from which we read the data
        # Get the variable
        try:
            self.var = obj.varMap[variable]
        except KeyError:
            error("ERROR: Variable: " + variable + " is not available for output for object named: " + obj.name + "\n" +
                  "       Valid variables are: " + str(list(obj.varMap.keys()))[:])

        self.expr = "self.obj." + self.var  # Expression to execute to evaluate the variable value

        self.file = open(self.location,
                         mode,
                         buffering=bufferSize)

        # Choose the writer
        value = eval(self.expr)
        if isinstance(value, np.floating) or isinstance(value, float) or isinstance(value, int):
            self.writer = NumWriter(obj, self.file)
        elif isinstance(value, dict):  # If list, split it)
            self.writer = DictWriter(obj, self.file, eval(self.expr))
        else:
            self.writer = ListWriter(obj, self.file)

    def write(self):
        value = eval(self.expr)
        self.writer.line(self.file, value)

    def close(self):
        self.file.close()

class Writer(ABC):
    """ Base class for writing data to files"""
    def __init__(self, obj, file):
        self.header(obj, file)

    def header(self, obj, file, data=None):
        try:
            name = obj.name
        except:
            name = "object"
        file.write('#' + name + '\n')

    @abstractmethod
    def line(self, file, value):
        pass

class ListWriter(Writer):
    def __init__(self, obj, file):
        super().__init__(obj, file)

    def line(self, file, value):
        file.write(" ".join(map(str, value)) + '\n')  # If list, split it

class NumWriter(Writer):
    def __init__(self, obj, file):
        super().__init__(obj, file)

    def line(self, file, value):
        file.write(str(value) + '\n')

class DictWriter:
    def __init__(self, obj, file, dict):
        self.header(obj, file, dict)

    def header(self, obj, file, dict):
        file.write('# ')
        for k, v in dict.items():
            file.write(k + 17 * ' ')
        file.write('\n')

    def line(self, file, dict):
        for k, v in dict.items():
            file.write(str(v) + 4 * " ")
        file.write('\n')


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

class ObjectRegistry:
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

