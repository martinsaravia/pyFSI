import pathlib
from pyFSI.post.anim import *


# Plot the eigenvalue progresssion
plt.figure()
xfilePath = pathlib.Path(__file__).parent.absolute()/'fsi'/'eigenValues.out'
yfilePath = pathlib.Path(__file__).parent.absolute()/'fsi'/'eigenValues.out'

plt.figure()
p = plotFromFile(xfilePath,
                 yfilePath,
                 xIndexes=np.s_[:, 0:10],
                 yIndexes=np.s_[:, 0:10],
				 imaginary=True)
p.axe.set_xlim([-100, 100])
p.axe.set_ylim([-600, 600])
p.axe.set_title("eigenLFBeam Test")
p.axe.set_xlabel("Eivenvalue real part")
p.axe.set_ylabel("Eigenvalue imaginary part")
plt.show()

