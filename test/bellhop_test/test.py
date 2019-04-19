import numpy as np
import sys
sys.path.append("../")
from os import system
from matplotlib import pyplot as plt
from pyat.objects import *
from pyat.readwrite import *



sd	=	20
#rd = [94.125, 99.755, 105.38, 111.00, 116.62, 122.25, 127.88, 139.12, 144.74, 150.38, 155.99, 161.62, 167.26, 172.88, 178.49, 184.12, 189.76, 195.38, 200.99, 206.62, 212.25]
rr	=	2.5

Z = np.arange(50, 51, 1)
X = np.arange(10, 11, 1)
print(X)

cw		=	1500
pw		=	1
aw		=	0
cb		=	1600
pb		=	1.5
ab		=	0.5

s = Source(sd)
r = Dom(X, Z)
pos = Pos(s,r)
pos.s.depth	= [sd]
pos.r.depth	 = Z
pos.r.range		=	X
pos.Nsd = 1
pos.Nrd = len(X)


depth = [0, 4000, 5000]

z1 = [0, 500, 1000, 1500, 2000, 3000, 4000]
# Layer 1
#z1		=	depth[:]	
alphaR	=	cw*np.array([1,.9, .8, .9, 1, 1, 1])
print(alphaR)
betaR	=	0.0*np.array([1]*len(z1))		
rho		=	pw*np.array([1]*len(z1))		
alphaI	=	aw*np.array([1]*len(z1))		
betaI	=	0.0*np.array([1]*len(z1))


# Layer 2
z2 = [4000, 5000]
#alphaR = 1600*[1, 1]
#betaR = [0, 0]
#rho

ssp1 = SSPraw(z1, alphaR, betaR, rho, alphaI, betaI)
print(ssp1.alphaR)


#	Sound-speed layer specifications
raw = [ssp1]
NMedia		=	1
Opt			=	'CVW'	
N			=	[0, 0]	
sigma		=	[0, 0]	
raw[0]
ssp = SSP(raw, depth, NMedia, Opt, N, sigma)

hs = HS()
Opt = 'A'
bottom = BotBndry(Opt, hs)
top = TopBndry('CVW')
bdy = Bndry(top, bottom)

plt.plot(ssp.sspf(np.linspace(0, 4000, 100)))
plt.show()

class Empty:
    def __init__(self):
        return

low = 1400
high = 1e9
cInt = CInt(low, high)
RMax = max(X)
freq = 50

# Beam params
run_type = 'R'
nbeams = 100
alpha = np.linspace(-5, 5, 100)
box = Box(1, 100)
deltas=20
beam = Beam(RunType=run_type, Nbeams=nbeams, alpha=alpha,Box=box,deltas=deltas)

write_env('py_env.env', 'BELLHOP', 'Pekeris profile', freq, ssp, bdy, pos, beam, cInt, RMax)
  

s = Source([sd])
ran =  np.arange(0,10, 10/1e3)
depth = np.arange(0,150,1)
r = Dom(ran, depth)

pos = Pos(s, r)

system("bellhop.exe py_env")
'''
system("/home/hunter/Downloads/at/bin/field.exe py_env")
[x,x,x,x,Pos1,pressure]= read_shd('py_env.shd')
#pressure = np.squeeze(pressure)
pressure = abs(pressure)
pressure = 10*np.log10(pressure / np.max(pressure))
print(np.shape(pressure))
levs = np.linspace(np.min(pressure), np.max(pressure), 20)
plt.contourf(Pos1.r.range, Pos1.r.depth,pressure[0, 0,:,:],levels=levs)
plt.gca().invert_yaxis()
plt.show()
'''
