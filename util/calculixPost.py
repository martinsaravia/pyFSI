import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pathlib
import os

FOAMRUN = pathlib.Path.home() / 'OpenFOAM/msaravia-7/run/openFoam7'

file = FOAMRUN / 'nagakura3.2/Solid/solid.dat'


data = np.loadtxt(file, comments=['displacements', 'velocities'])
dt = 1E-3

disp = data[2:-1:2, 2]
vel = data[1:-1:2, 2]

time = np.arange(0, len(disp)*dt, dt)
axe = plt.gca()
axe.set_title(r'Time History', fontsize=14)
axe.set_xlabel('Time (s)', fontsize=14)
axe.set_ylabel('Displacement (m)',fontsize=14)
axe.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
plt.rc('text', usetex=True)
plt.rc('font', **{
    'family': 'serif',
    'serif': ['Times']
})
#print(len(disp), len(vel))

# plt.scatter(disp, vel)
plt.plot(time, disp)
plt.show()
