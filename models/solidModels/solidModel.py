from abc import ABCMeta, abstractmethod
import os, pathlib

class solidModel(metaclass=ABCMeta):
    def __repr__(self):
        return 'solidModel Abstract Class'

    def __init__(self, control, mesh, dict, name='NN'):
        self._control = control
        self._name = name
        self._mesh = mesh
        self._dict = dict

    def name(self):
        return self._name

    def mesh(self):
        return self._mesh

    def dict(self):
        return self._dict
