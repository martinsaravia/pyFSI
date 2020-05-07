
from solvers.caseControl import *
import matplotlib.pyplot as plt, numpy as np
import matplotlib.cm as cm
solution = runCase("caltechHarvester01", parametric=True)


# Profiled run
# import cProfile, pstats
# cProfile.run('runCase("caltechHarvester00")', 'profileRes')
# prof = pstats.Stats('profileRes')
# prof.sort_stats('tottime').print_stats(10)
# prof.sort_stats('cumtime').print_stats(10)
modes = 2
esize = modes*2+2
colors = cm.rainbow(np.linspace(0, 1, esize))
for key, val in solution.items():
    esystems = [s.eigen() for s in solution[key]]

    for es in esystems:
        plt.scatter(np.real(es.evalues()), np.imag(es.evalues()), color=colors)
    # plt.xlim(-100,100)
    # plt.ylim(0,600)
    plt.show()
