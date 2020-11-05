#!/usr/bin/env python3
import math
import sys

import pycalculix as pyc

# We'll be modeling a masonry gravity dam, the Beetaloo dam in Australia
# Problem constants
proj_name = 'beam'
grav = 9.81 # m/s^2


# set whether or not to show gui plots
show_gui = True
if '-nogui' in sys.argv:
    show_gui = False
# set element shape
eshape = 'quad'
if '-tri' in sys.argv:
    eshape = 'tri'


# make model
model = pyc.FeaModel(proj_name)
# this sets dist units to meters, labels our consistent units
model.set_units('m')

# make part, coordinates are x, y = radial, axial
part = pyc.Part(model)
lines = []
part.goto(0.0,0.0)
[L1,p1,p2] = part.draw_line_to(10.0,0.0)
lines.append(L1)
[L1,p1,p2] = part.draw_line_to(10.0,-0.7)
lines.append(L1)
[L1,p1,p2] = part.draw_line_to(0.0,-0.7)
lines.append(L1)
part.draw_line_to(0, 0)
model.plot_geometry(proj_name+'_geom', display=show_gui)

# set part material
mat = pyc.Material('concrete')
mat.set_mech_props(2300, 30000*(10**6), 0.2)
model.set_matl(mat, part)

# set the element type, line division, and mesh the database
thickness = 1.0
eshape = 'quad'
model.set_eshape(eshape, 2)
model.set_etype('plstrain', part, thickness)
model.get_item('L0').set_ediv(50)
model.get_item('L1').set_ediv(4)
model.get_item('L2').set_ediv(50)
model.get_item('L3').set_ediv(4)
print(lines)
# mesh with 1.0 or less fineness, smaller is finer
model.mesh(size=0.2, meshmode='fineness', mesher='gmsh')
# plot the part elements
model.plot_elements(proj_name+'_elem', display=show_gui)

# set loads and constraints
model.set_constr('fix', part.bottom, 'x')
model.set_constr('fix', part.bottom, 'y')
model.plot_constraints(proj_name+'_constr', display=show_gui)

# make model and solve it
prob = pyc.Problem(model, 'struct')
prob.solve()

# query results and store them
disp = False  # turn off display plotting
# store the fields to write
fields = 'Seqv,Sx,Sy,Sz,S1,S2,S3,ux,uy,utot'
fields = fields.split(',')
for field in fields:
    fname = proj_name+'_'+field
    prob.rfile.nplot(field, fname, display=disp)
smax = prob.rfile.get_nmax('Seqv')
[fx, fy, fz] = prob.rfile.get_fsum(model.get_item('L9'))
print('Seqv_max= %3.2f' % (smax))
print('Reaction forces (fx,fy,fz) = (%12.10f, %12.10f, %12.10f)' % (fx, fy, fz))
