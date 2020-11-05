from post.plots import *
from post.anim import *
import importlib, json, pathlib, sys

fileName = "beamNoFlow.out"

# filePath = pathlib.Path(__file__).parent.absolute()/fileName
#
# pFile(filePath, 0)

# filePath = pathlib.Path(__file__).parent.absolute()/'solid'/'yMid.out'
# p = pFile(filePath,  indexes=[-1], time='all')
# #
# plt.figure()
# filePath = pathlib.Path(__file__).parent.absolute()/'fluid'/'Q0.out'
# p = pFile(filePath,  indexes=[0, 1], time='all')

# plt.figure()
# filePath = pathlib.Path(__file__).parent.absolute()/'solid'/'yMid.out'
# p = pFile(filePath,  indexes=[-1], time='all')
#
# plt.figure()
# filePath = pathlib.Path(__file__).parent.absolute()/'solid'/'dyMid.out'
# p = pFile(filePath, indexes=[-1], time='all')
#
# plt.figure()
# filePath = pathlib.Path(__file__).parent.absolute()/'solid'/'ddyMid.out'
# p = pFile(filePath, indexes=[-1], time='all')
#
# plt.figure()
# filePath = pathlib.Path(__file__).parent.absolute()/'fluid'/'Q0.out'
# p = pFile(filePath,
#           indexes=[0, 1],
#           time=[0, 3100],
#           factor=(1.0/5E-3))
#
# p.axe.set_xlim([0, 3.3333])
# p.axe.set_ylim([0, 5])
# p.axe.set_title("Nagakura Test")
# p.axe.set_xlabel("Pressure difference (Pa)")
# p.axe.set_ylabel("Inlet velocity (m/s)")
# p.axe.legend(['Top channel', 'Bottom channel'])
# # plt.show()
# anim = pAnimation(p, name='response', interval=10, frames=3100)
# anim.save('FlowRate', fps=120 )

plt.figure()
filePath = pathlib.Path(__file__).parent.absolute()/'solid'/'yMid.out'
p = pShape(filePath, time=3100)
p.axe.set_title("Nagakura Test")
p.axe.set_xlabel("x (m)")
p.axe.set_ylabel("Displacement (mm)")
p.axe.set_xlim([0, 0.22])
anim = pAnimation(p, name='response', interval=10, frames=3100)
anim.save('Nagakura2', fps=120 )


