import pathlib

from post.anim import *


# plt.figure()
# xfilePath = pathlib.Path(__file__).parent.absolute()/'fluid'/'pIn.out'
# yfilePath = pathlib.Path(__file__).parent.absolute()/'fluid'/'v0.out'
# p = pFile(xfilePath,      yfilePath,
#           xIndexes='all', yIndexes=[0, 1])
# p.axe.set_xlim([0, 5000])
# p.axe.set_ylim([0, 300])
# p.axe.set_title("Vaughan SN 01")
# p.axe.set_xlabel("Pressure difference (Pa)")
# p.axe.set_ylabel("Inlet velocity (m/s)")
# p.axe.legend(['Top channel', 'Bottom channel'])
# plt.show()
#
#
# plt.figure()
# xfilePath = pathlib.Path(__file__).parent.absolute()/'fluid'/'time.out'
# yfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'yMid.out'
# p = pFile(xfilePath,      yfilePath,
#           xIndexes='all', yIndexes=[-1])
# p.axe.set_xlim([0, 15])
# p.axe.set_ylim([-1E-2, 1E-2])
# p.axe.set_title("Vaughan SN 01")
# p.axe.set_xlabel("Time(s)")
# p.axe.set_ylabel("Displacement (m)")
# p.axe.legend(['Top channel', 'Bottom channel'])
# plt.show()
#
#
# plt.figure()
# xfilePath = pathlib.Path(__file__).parent.absolute()/'fluid'/'time.out'
# yfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'dyMid.out'
# p = pFile(xfilePath,      yfilePath,
#           xIndexes='all', yIndexes=[-1])
# p.axe.set_xlim([0, 15])
# p.axe.set_ylim([-5, 5])
# p.axe.set_title("Vaughan SN 01")
# p.axe.set_xlabel("Time (s)")
# p.axe.set_ylabel("Tip velocity (m/s)")
# p.axe.legend(['Top channel', 'Bottom channel'])
# plt.show()
#
# plt.figure()
# tfilePath = pathlib.Path(__file__).parent.absolute()/'fluid'/'time.out'
# xfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'x.out'
# yfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'yMid.out'
# p = pShape(xfilePath, yfilePath, time=3100)
# p.axe.set_title("Vaughan 01")
# p.axe.set_xlabel("x (m)")
# p.axe.set_ylabel("Displacement (m)")
# p.axe.set_xlim([0, 0.12])
# p.axe.set_ylim([-0.01, 0.01])
#
# anim = pAnimation(p, tfilePath, name='response', interval=500)
# anim.save('Vaughan01', fps=120 )


plt.figure()
filePath = pathlib.Path(__file__).parent.absolute()/'magnetic'/'Vaughan01_5mm.flux'
p = pFlux(filePath, deriv=True)
# p.axe.set_xlim([0, 15])
# p.axe.set_ylim([-5, 5])
p.axe.set_title("Vaughan SN 01")
p.axe.set_xlabel("x (m)")
p.axe.set_ylabel("Magnetic Flux Derivative")
plt.show()
