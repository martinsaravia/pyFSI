# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt, numpy as np
import matplotlib.cm as cm
import scipy.signal as sps
from matplotlib.animation import FuncAnimation
import pandas as pd
from multiprocessing import Pool

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "figure.figsize": [16, 9],
    "savefig.format": "png",
    "animation.writer": "ffmpeg",
    "animation.ffmpeg_args": ['-c:v', 'libx264', '-crf', '0'],
    "animation.bitrate": 24000,
    "animation.frame_format": "png",
    "legend.loc": 'upper left'
})


# plt.tight_layout()


def loadFiles(fileList):
    n = len(fileList)
    pool = Pool(2)
    result = pool.map(np.loadtxt, fileList)
    pool.close()
    return result

    # self.axe.set_xticks([])
    # self.axe.set_yticks([])


# Load files in parallel
class fsiPlot():
    def __init__(self):
        self.fig = plt.gcf()
        self.axe = plt.gca()
        self.xData = None
        self.yData = None
        plt.grid(b=True)

        # plt.rcParams.update({
        #     "text.usetex": True,
        #     "font.family": "sans-serif",
        #     "font.sans-serif": ["Helvetica"]})
        ## for Palatino and other serif fonts use:

    def set(self,
            fig=False,
            axe=False,
            title=".",
            xlim=False,
            ylim=False,
            xLabel=".",
            yLabel="."):

        if fig:
            self.fig = fig
        if axe:
            self.axe = axe
        if xlim:
            self.axe.set_xlim(xlim[0], xlim[1])
        if ylim:
            self.axe.set_ylim(ylim[0], ylim[1])

        # Set the labels
        self.axe.set_title(title)
        self.axe.set_xlabel(xLabel)
        self.axe.set_ylabel(yLabel)


class plotFromFile(fsiPlot):
    def __init__(self, xFile, yFile, xIndexes=None, yIndexes=None):
        super().__init__()
        data = loadFiles([xFile, yFile])

        self.P = []

        # Number of plots
        try:
            nPlots = abs(yIndexes[1].stop - yIndexes[1].start)
        except:
            nPlots = 1

        # Number of x Axes
        try:
            xAxes = abs(xIndexes[1].stop - xIndexes[1].start)
        except:
            xAxes = 1

        # Solve the problem of different lengths for an incomplete run
        # if xAxes != 1:
        #     if xIndexes[0].stop == -1:
        #         if xIndexes[0].stop > yIndexes[0].stop:
        #             xIndexes[0].stop = yIndexes[0].stop
        #         if xIndexes[0].stop < yIndexes[0].stop:
        #             yIndexes[0].stop = xIndexes[0].stop
        # else:
        #     print(xIndexes)
        #     if xIndexes.stop == -1:
        #         if xIndexes.stop > yIndexes[0].stop:
        #             xIndexes.stop = yIndexes[0].stop
        #         if xIndexes.stop < yIndexes[0].stop:
        #             yIndexes[0].stop = xIndexes.stop

        self.xData = np.array(data[0])[xIndexes]
        self.yData = np.array(data[1])[yIndexes]

        if xAxes == 1 and nPlots == 1:
            curve, = self.axe.plot(self.xData, self.yData)
            self.P.append(curve)

        else:
            if xAxes == 1:
                for i in range(nPlots):
                    curve, = self.axe.plot(self.xData, self.yData[:, i])
                    self.P.append(curve)
            else:
                for i in range(nPlots):
                    curve, = self.axe.plot(self.xData[:, i], self.yData[:, i])
                    self.P.append(curve)

    def update(self, i):
        for c in self.yIndexes:
            self.P[c].set_xdata(self.xData[xIndexes[0]:i])
            self.P[c].set_ydata(self.yData[xIndexes[0]:i, c])
        return self.P


class pShape(fsiPlot):
    def __init__(self, xfile, yfile, time):
        super().__init__()
        self.time = np.loadtxt(xfile)
        self.data = np.loadtxt(yfile)
        self.P = self.axe.plot(self.time, self.data[time, :], linewidth=3)

    def update(self, i):
        self.P[0].set_ydata(self.data[i, :])
        return self.P[0]


class pScatter(fsiPlot):
    def __init__(self, solution):
        # pidx is the parameter index of the solution list
        # tidx is the time  index
        super().__init__()
        if solution[0][0].execution()['isParametric']:
            para = solution[0][0].execution()['parameters']
        else:
            para = {'p': [1], 'pi': [1], 'steps': 1}

        time = solution[0][0].execution()['time']
        esize = solution[0][0].ES.size()
        self.colors = cm.rainbow(np.linspace(0, 1, esize))
        self.data = np.zeros((para['steps'], time['steps'], esize, 2))
        self.evec = np.zeros((para['steps'], time['steps'], esize, esize),
                             dtype=complex)
        self.P = self.axe.scatter([], [], edgecolors='k', facecolors='r', s=10)
        for pi, ppf in enumerate(para['p']):  # Parameter loop
            for ti, t in enumerate(time['t']):  # Time loop
                evalues = solution[pi][ti].ES.evalues()
                evectors = solution[pi][ti].ES.evectors()
                # Store real and imaginary parts of eigenvalues
                self.data[pi, ti, :, 0] = np.real(evalues)
                self.data[pi, ti, :, 1] = np.imag(evalues)
                self.evec[pi, ti, :, :] = evectors

        self.esize = esize
        self.tsize = len(time['t'])

    def plot(self, i):
        x = self.data[i, :, :, 0].flatten()
        y = self.data[i, :, :, 1].flatten()
        self.P = self.axe.scatter(x, y, edgecolors='k', facecolors='r', s=5)
        self.fig.show()

    def update(self, i):
        x = self.data[i, :, :, 0].flatten()
        y = self.data[i, :, :, 1].flatten()
        data = np.array([x, y]).T  # The transpose MUST be called
        self.P.set_offsets(data)
        self.fig.show()


# Class for plotting magnetic flux function taken from FEMM with lua script
class pFlux(fsiPlot):
    def __init__(self, filePath, deriv=True):
        super().__init__()
        self.data = np.genfromtxt(filePath, delimiter=",", skip_header=3)

        # Flux
        self.aflx = np.zeros((len(self.data[:, 0]), 3))
        self.aflx[:, 0] = self.data[:, 0]
        self.aflx[:, 1] = sps.savgol_filter(self.data[:, 1], 11, 2, deriv=0)
        self.aflx[:, 2] = self.data[:, 1]
        # Derivative of the flux
        self.dflx = np.zeros((len(self.aflx[:, 0]), 2))  # Derivative of average flux
        self.dflx[:, 0] = self.aflx[:, 0]  # LENGTH COORDINATE
        temp = sps.savgol_filter(self.data[:, 1], 101, 2, deriv=1)
        # Necesito el dx para dividir la derivada porqeu sale calculada tomando dx=1 por defecto
        dx = self.aflx[1, 0] - self.aflx[0, 0]
        # Ojo que tira la derivada cambiada de signo
        self.dflx[:, 1] = -sps.savgol_filter(temp, 51, 2, deriv=0) / dx

        # Choose the flux function or the derivative
        if deriv:
            curve, = self.axe.plot(self.dflx[:, 0], self.dflx[:, 1])
        else:
            curve, = self.axe.plot(self.aflx[:, 0], self.aflx[:, 1])

        self.P = []

        self.P.append(curve)


class pBifurcation(fsiPlot):
    def __init__(self, solution, mode0, mode1, var='freq'):
        super().__init__()
        self.P = self.axe.scatter([], [], edgecolors='k', facecolors='r', s=10)
        self.bif = []  # Bifurcation points list
        size = len(solution[0][0].ES.evalues())

        for sol in solution:  # Loop through the parameters
            oldReals = -np.ones(size)  # Assume stability
            searchModes = range(mode0, mode1)
            for fsi in sol:
                newReals = np.real(fsi.ES.evalues())
                newImags = np.imag(fsi.ES.evalues())
                signs = np.sign(newReals) + np.sign(oldReals)

                # Extract all unstable states
                for mode in searchModes:
                    if newReals[mode] > 0 and newImags[mode] > 0:
                        xVar = fsi.execution()['parameters']['pi']
                        if var == 'frequency':
                            yVar = newImags[mode]
                        if var == 'flowRate':
                            yVar = fsi.flow().qx0
                        elif var == 'velocity':
                            yVar = fsi.flow().vRef
                        else:
                            yVar = None
                            print("ERROR: No yVar defined...")
                        self.bif.append(bifurcationPoint(xVar, yVar, mode))

                # Extrac the stability boundary
                # for mode in searchModes:
                #     if signs[mode] == 0 and np.sign(oldReals[mode]) != 0 and newImags[mode] > 0:
                #         self.bif.append(bifurcationPoint(t.pi, newImags[mode], mode))
                #         #searchModes = [mode] # If I found one unstable, follow this
                oldReals = newReals

        # Get the data from the bifurcation point
        self.xAxis = []
        self.yAxis = []
        for p in self.bif:
            self.xAxis.append(p.xParameter)
            self.yAxis.append(p.yParameter)

    def plot(self):
        self.axe.scatter(self.xAxis,
                         self.yAxis,
                         edgecolors='k',
                         facecolors='r',
                         s=10)
        plt.show()


class bifurcationPoint():
    def __init__(self, xPar, yPar, neval):
        self.xParameter = xPar
        self.yParameter = yPar
        self.modeNumber = neval


