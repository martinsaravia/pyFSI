# Run the case
from execution import execution
execution.cleanCase()
case = execution.runCase("eigenLFBeam")


# Plots
import pathlib
from post.anim import *
# Load the result files
xfilePath = pathlib.Path(__file__).parent.absolute()/'fsi'/'realValues.out'
yfilePath = pathlib.Path(__file__).parent.absolute()/'fsi'/'imagValues.out'

plt.figure()
p = plotFromFile(xfilePath,
                 yfilePath,
                 xIndexes=np.s_[:, 0:10],
                 yIndexes=np.s_[:, 0:10])
p.axe.set_xlim([-100, 100])
p.axe.set_ylim([-600, 600])