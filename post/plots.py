# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt, numpy as np
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation

plt.rcParams.update({
    "font.family": "serif",
    "font.size"   : 12,
    "figure.figsize": [16, 9],
    "savefig.format": "png",
    "animation.writer": "ffmpeg",
    "animation.ffmpeg_args": ['-c:v', 'libx264', '-crf', '0'],
    "animation.bitrate": 24000,
    "animation.frame_format": "png",
    "legend.loc": 'upper left'
})

# plt.tight_layout()



class fsiPlot():
    def __init__(self):
        self.fig = plt.gcf()
        self.axe = plt.gca()
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


        # self.axe.set_xticks([])
        # self.axe.set_yticks([])

class pFile(fsiPlot):
    def __init__(self, file, indexes=None, time='all', factor=1):
        super().__init__()
        self.factor = factor
        self.indexes = indexes
        self.time = time
        self.data = np.loadtxt(file)
        self.P = []
        self.x = np.linspace(0, 3.3333, time[1])
        if time == 'all':
            for i in indexes:
                self.P = self.axe.plot(factor * self.data[:, i])
        else:
            for i in indexes:
                curve, = self.axe.plot(self.x, factor * self.data[time[0]:time[1], i])
                self.P.append(curve)

    def update(self, i):
        for c in self.indexes:
            self.P[c].set_xdata(self.x[self.time[0]:i])
            self.P[c].set_ydata(self.factor * self.data[self.time[0]:i, c])
        return self.P

class pShape(fsiPlot):
    def __init__(self, file, time):
        super().__init__()
        self.time = time
        self.data = np.loadtxt(file)
        x = np.linspace(0, 0.2, 21)
        self.P = self.axe.plot(x, self.data[time, :], linewidth=3)

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
        self.evec = np.zeros((para['steps'], time['steps'], esize, esize), dtype=complex)
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









#         

#     if not update:
#         fig = plt.figure()
#         ax = plt.axes()
#         colors = cm.rainbow(np.linspace(0, 1, esize))
#         plt.xlim(-100,100)
#         plt.ylim(-750,750)
        
#     # Gather the eigensystems
#     for key, val in solution.items():
#         esystems = [s.ES for s in solution[key]]
        
        

#     for es in esystems:
#         plt.scatter(np.real(es.evalues()), np.imag(es.evalues()), color=colors)
#         print(np.real(es.evalues()[0:esize]))
#     print("---------------------------------")
#     # plt.xlim(-100,100)
#     # plt.ylim(0,600)
#     plt.show()
    
    

#     modes = 2
# esize = modes*2+2
# colors = cm.rainbow(np.linspace(0, 1, esize))




# for key, val in solution.items():
#     esystems = [s.ES for s in solution[key]]

#     for es in esystems:
#         plt.scatter(np.real(es.evalues()), np.imag(es.evalues()), color=colors)
#         print(np.real(es.evalues()[0:esize]))
#     print("---------------------------------")
#     # plt.xlim(-100,100)
#     # plt.ylim(0,600)
#     plt.show()