import numpy as np
from matplotlib import pyplot as plt
from env.env.envs import factory
from pyat.pyat.env import Beam, Box
from pyat.pyat.readwrite import read_arrivals_asc, read_shd

'''
Description:

Author: Hunter Akins
'''

builder = factory.create('deepwater')
dw_env= builder()

freq = 1400
zs = 500
zr = np.linspace(0, 5000, 20)
dz = 20
zmax = 5000
dr = 1
rmax = 100*1e3

dw_env.add_source_params(freq, zs, zr)
dw_env.add_field_params(dz, zmax, dr, rmax)




folder = 'at_files/'
fname = 'dw'


"""
Impulse response
"""

run_type='A'
nbeams=131
alpha=np.linspace(-20, 20, nbeams)
box=Box(zmax+100, rmax*1e-3)
deltas = 12
beam = Beam(RunType=run_type, Nbeams=nbeams, alpha=alpha,box=box,deltas=deltas)
dw_env.run_model('bellhop', folder, fname, beam=beam)
arrivals, pos = read_arrivals_asc(folder+fname)


"""
Shd plot
"""
beam.RunType = 'I'
dw_env.run_model('bellhop', folder, fname, beam=beam, zr_range_flag=False, zr_flag=False)

[x,x,x,x,ppos, p] = read_shd(folder + 'dw.shd')
print(p.shape)
p = abs(p)
p = 10*np.log10(p/np.max(p))
print(p)
levs = np.linspace(-60, 0, 20)
plt.contourf(np.squeeze(p), levels=levs)
plt.gca().invert_yaxis()
plt.show()
#print(len(arrivals))
#r_arrivals = np.array(arrivals[1])
#print(r_arrivals)
#arr_times = r_arrivals[:,1]
#print(arr_times)
#
#def impulse_response(t):
#    if min([abs(t-x) for x in arr_times]) < 1e-3:
#        return 1
#    else:
#        return 0 
#t = np.linspace(3, 3.4, 1000)
#plt.plot(t, np.array([impulse_response(x) for x in t]))
#plt.show()


