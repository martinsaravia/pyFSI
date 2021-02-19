from abc import ABCMeta, abstractmethod

class fluidModel(metaclass=ABCMeta):
    def __repr__(self):
        return 'fluidModel Abstract Class'

    def __init__(self, fluid):
        self._rho = fluid['rho']
        self._nu = fluid['nu']

    @abstractmethod
    def getP(self, loc):
        pass
