import pathlib

from pyFSI.post.anim import *

# Pressure vs velocity
plt.figure()
# Set the path to the files
xfilePath = pathlib.Path(__file__).parent.absolute()/'flow'/'pIn.out'
yfilePath = pathlib.Path(__file__).parent.absolute()/'flow'/'v0.out'
p = plotFromFile(xfilePath,
                 yfilePath,
                 xIndexes=np.s_[:],
                 yIndexes=np.s_[:, 0:2])
p.axe.set_xlim([0, 6000])
p.axe.set_ylim([0, 150])
p.axe.set_title("Vaughan SN 01")
p.axe.set_xlabel("Pressure difference (Pa)")
p.axe.set_ylabel("Inlet velocity (m/s)")
p.axe.legend(['Top channel', 'Bottom channel'])
plt.show()
#
#Time vs displacement plot
plt.figure()
xfilePath = pathlib.Path(__file__).parent.absolute()/'flow'/'time.out'
yfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'yMid.out'
p = plotFromFile(xfilePath,
                 yfilePath,
                 xIndexes=np.s_[0:28000],
                 yIndexes=np.s_[0:28000, -1])
p.axe.set_xlim([0, 5])
p.axe.set_ylim([-1E-2, 1E-2])
p.axe.set_title("Vaughan SN 01")
p.axe.set_xlabel("Time(s)")
p.axe.set_ylabel("Displacement (m)")
p.axe.legend(['Top channel', 'Bottom channel'])
plt.show()
#
#
#
## Beam deformed shape animation
plt.figure()
tfilePath = pathlib.Path(__file__).parent.absolute()/'flow'/'time.out'
xfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'x.out'
yfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'yMid.out'
p = pShape(xfilePath, yfilePath, time=3100)
p.axe.set_title("Vaughan 01")
p.axe.set_xlabel("x (m)")
p.axe.set_ylabel("Displacement (m)")
p.axe.set_xlim([0, 0.12])
p.axe.set_ylim([-0.01, 0.01])
anim = pAnimation(p, tfilePath, name='response', interval=2000)
anim.save('transientLFBeam', fps=120 )
