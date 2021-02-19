

file = openfile ('Vaughan01_5mm.flux', 'w')
ri = 0.45
re = 0.55
z0 = -0.2
z1 =  0.2
dz = z1 - z0
samp = 1000
step = dz / samp
z = z0
write(file, '** -------------------------------------------------------------------------', '\n' )
write(file, '**   re = ', re, '   z0 = ', z0, '   dz = ', dz, '   samp = ', samp, '   step = ', step, '\n' )
write(file, '** -------------------------------------------------------------------------', step, '\n' )
for i = 1, samp, 1 do
	mo_addcontour(ri, z)	
	mo_addcontour(re, z)
	inte, favg = mo_lineintegral(0)
	print(i,z,favg)
	mo_clearcontour()
	write(file, z, ',', favg, '\n' )
	z = z + step
end 
closefile(file)

