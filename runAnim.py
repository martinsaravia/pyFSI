
from solvers.caseControl import runCase
import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from post.plots import pScatter, pBifurcation
from post.prints import *
import importlib

import pandas as pd

solution = runCase("caltechHarvester02")

import post
from post.plots import pScatter, pBifurcation
importlib.reload(post.plots)
p0 = pBifurcation(solution)
p0.plot()


#
plt.close('all')
p = pScatter(solution)
for i in range(len(solution[0])):


printStateMatrix(solution[0][0])


















# outfile = open('temp', 'wb')
# pickle.dump(solution, outfile)
# outfile.close()
# outfile = open('temp', 'rb')
# sol2 = pickle.load(outfile)

# print(sol2)
# Profiled run
# import cProfile, pstats
# cProfile.run('runCase("caltechHarvester00")', 'profileRes')
# prof = pstats.Stats('profileRes')
# prof.sort_stats('tottime').print_stats(10)
# prof.sort_stats('cumtime').print_stats(10)
# frames = len(solution.items())
# modes = 2
# esize = modes*2+2
# #
# # fig = plt.figure(figsize=(7, 7))
# # ax = fig.add_axes([0, 0, 1, 1], frameon=False)


# fig = plt.figure()
# ax = plt.axes()
# colors = cm.rainbow(np.linspace(0, 1, esize))
# plt.xlim(-100,100)
# plt.ylim(-750,750)

# steps = len(solution.values())
# nesys = len(list(solution.values())[0])
# D = np.zeros((steps, nesys, esize, 2))
# colors = cm.rainbow(np.linspace(0, 1, esize))

# mode0 = 0
# mode1 = nesys
# i = 0
# for key, val in solution.items():
#     esystems = [s.eigen() for s in solution[key]]
#     for idx, es in enumerate(esystems):
#         D[i, idx, :, 0] = np.real(es.evalues())[:]
#         D[i, idx, :, 1] = np.imag(es.evalues())[:]
#     i += 1

# x = []
# y = []
# for i in range(nesys):
#     x.append(D[0, i, :, 0])
#     y.append(D[0, i, :, 1])
# xflat = [item for sublist in x for item in sublist]
# yflat = [item for sublist in y for item in sublist]

# scat = ax.scatter(xflat,yflat, edgecolors='k', facecolors='r')
# data = np.zeros((len(xflat),2))

# def update(frame):
#     x = []
#     y = []
#     for i in range(nesys):
#         x.append(D[frame, i, :, 0].tolist())
#         y.append(D[frame, i, :, 1].tolist())
#     data[:,0] = [item for sublist in x for item in sublist]
#     data[:,1] = [item for sublist in y for item in sublist]
#     scat.set_offsets(data)


# ani = FuncAnimation(fig, update, interval=200, frames=100)
# plt.show()



# plt.show()
