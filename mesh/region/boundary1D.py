from abc import ABCMeta, abstractmethod
import numpy as np

class boundary1D(metaclass=ABCMeta):

    def __repr__():
        return 'boundary1D'

    def __init__(self, mesh, name):
        self._name = name
        self._mesh = mesh

        # First derivative
        edgeOrder = 2
        self._dy = np.gradient(self.y(), self._mesh.x(), edge_order=edgeOrder)

        # Powers
        self._y3 = self.y()**3
        self._y4 = self.y()**4


    # Abstract methods (Pure virtual function, must be implemented in subclasses)
    @abstractmethod
    def update(self):
        pass

    @abstractmethod
    def isFlexible(self):
        pass

    @abstractmethod
    def y(self):
        pass

    @abstractmethod
    def ytop(self):
        pass

    @abstractmethod
    def ybot(self):
        pass

    # Getters
    def __call__(self):
        return self.y()

    def name(self):
        return self._name

    def mesh(self, name):
        return self._mesh

    def dy(self):
        return self._dy

    def y3(self):
        return self._y3

    def y4(self):
        return self._y4
