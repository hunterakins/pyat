import numpy as np
import sys
sys.path.append("../")
from os import system
from matplotlib import pyplot as plt
from pyat.pyat.env import *
from pyat.pyat.readwrite import *




sd	=	20
rd = [94.125, 99.755, 105.38, 111.00, 116.62, 122.25, 127.88, 139.12, 144.74, 150.38, 155.99, 161.62, 167.26, 172.88, 178.49, 184.12, 189.76, 195.38, 200.99, 206.62, 212.25]
rr	=	2.5

Z = np.arange(0, 100, .5)
X = np.arange(0, 10, .01)

cw		=	1500
pw		=	1
aw		=	0
cb		=	1600
pb		=	1.8
ab		=	0.2

s = Source(sd)
r = Dom(X, Z)
pos = Pos(s,r)
pos.s.depth	= [sd]
pos.r.depth	 = Z
pos.r.range		=	X
pos.Nsd = 1
pos.Nrd = len(rd)



bottom_depth = 100
depth = [0, bottom_depth] 
# Layer 1
z1		=	depth[0:2]	
z1 = np.linspace(depth[0], bottom_depth, 1000)
alphaR	=	cw*np.ones(z1.shape)
betaR	=	0.0*np.ones(z1.shape)
rho		=	pw*np.ones(z1.shape)
alphaI	=	aw*np.ones(z1.shape)
betaI	=	0.0*np.ones(z1.shape)

ssp1 = SSPraw(z1, alphaR, betaR, rho, alphaI, betaI)

#	Sound-speed layer specifications
raw = [ssp1]
NMedia		=	1
Opt			=	'CVW'	
N			=	[z1.size]
sigma		=	[.5,.5]	 # roughness at each layer. only effects attenuation (imag part)
ssp = SSP(raw, depth, NMedia, Opt, N, sigma)


hs = HS(alphaR=cb, betaR=0, rho = pb, alphaI=ab, betaI=0)
Opt = 'A~'
bottom = BotBndry(Opt, hs)
top = TopBndry('CVW')
bdy = Bndry(top, bottom)


class Empty:
    def __init__(self):
        return

cInt = Empty()
cInt.High = cb
cInt.Low = 0 # compute automatically
RMax = max(X)
freq = 100
write_env('py_env.env', 'KRAKEN', 'Pekeris profile', freq, ssp, bdy, pos, [], cInt, RMax)
  

s = Source([sd])
ran =  np.arange(0,10, 10/1e3)
depth = np.arange(0,1.5*bottom_depth,1)
r = Dom(ran, depth)

pos = Pos(s, r)

write_fieldflp('py_env', 'R', pos)
system("krakenc.exe py_env")
fname = 'py_env.mod'
options = {'fname':fname, 'freq':0}
modes = read_modes(**options)
print(modes)
print(modes.k)
delta_k = np.max(modes.k.real) - np.min(modes.k.real)
print('range cell size', 2*np.pi/delta_k)
bandwidth = delta_k * 2.5 / 2 / (2*np.pi)
print('bandwidth', bandwidth )
print('coh time', 1/bandwidth)
print('cell cross time', 2*np.pi / delta_k / 2.5)
#figs = modes.plot()
#plt.show()
system("field.exe py_env")
[x,x,x,x,Pos1,pressure]= read_shd('py_env.shd')
#pressure = np.squeeze(pressure)
pressure = abs(pressure)
pressure = 10*np.log10(pressure / np.max(pressure))
print(np.shape(pressure))
levs = np.linspace(np.min(pressure), np.max(pressure), 20)
plt.contourf(Pos1.r.range, Pos1.r.depth,(pressure[0, 0,:,:]),levels=levs)
plt.gca().invert_yaxis()
plt.show()

