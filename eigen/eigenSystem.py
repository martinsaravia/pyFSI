import numpy as np

class eigenSystem(object):
    def __init__(self):
        pass

    def form(self, eigenValueList, eigenVectorList):
        self.__system =  []
        for i in eigenValueList:
            self.__system.append( [eigenValueList[i], eigenValueList[i]] )

    def eVector(self, i):
        return self.__system[i,1]

    def eValue(self, i):
        return self.__system[i,0]
