
from solvers.caseControlParallel import runCase
import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib.cm as cm
from post.plots import pScatter, pBifurcation
from post.anim import pAnimation
from post.prints import *
import importlib
import pandas as pd

#solution = runCase("TosiFig3.1")
solution = runCase("Nagakura")
# solution = runCase("TosiFig3.1")

sol0 = solution[1][0]
import post
from post.plots import pScatter, pBifurcation
import post
from post.plots import pScatter, pBifurcation
importlib.reload(post.plots)
p0 = pBifurcation(solution, 0, 10)
p0.set(xlim=[1, 4], ylim=[0, 5E-2])
p0.plot()


import post
from post.plots import pScatter
import post
from post.plots import pScatter
p = pScatter(solution)
p.set(yLabel='Flow Rate - $\mathbf{Q}$',
      xLabel='Throat size - $\mathit{h}$',
      title=r'FSI Harvester - Eigenvalue progression',
      xlim=[-100, 100],
      ylim=[-750, 750])
p.plot(1)

# anim = pAnimation(p, name="eigenSolution", frames=len(solution))

# # anim.save('TosiFig32')



# importlib.reload(post.prints)
printNumbers(solution[0][-1])
# #
# #
# # plt.close('all')
# for i, sol in enumerate(solution):
#     p = pScatter(solution)
#     p.plot(i)

# printStateMatrix(solution[0][0])















