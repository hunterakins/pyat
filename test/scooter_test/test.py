import numpy as np
import sys
sys.path.append("../")
from os import system
from matplotlib import pyplot as plt
from pyat.pyat.env import *
from pyat.pyat.readwrite import *



sd	=	30
rd = [94.125, 99.755, 105.38, 111.00, 116.62, 122.25, 127.88, 139.12, 144.74, 150.38, 155.99, 161.62, 167.26, 172.88, 178.49, 184.12, 189.76, 195.38, 200.99, 206.62, 212.25]
rr	=	2.5

Z = np.arange(.5, 150.5, .5)
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



depth = [0, 100, 150]
# Layer 1
z1		=	depth[0:2]	
alphaR	=	cw*np.array([1,1])
betaR	=	0.0*np.array([1,1])		
rho		=	pw*np.array([1,1])		
alphaI	=	aw*np.array([1,1])		
betaI	=	0.0*np.array([1,1])

ssp1 = SSPraw(z1, alphaR, betaR, rho, alphaI, betaI)
arr = np.linspace(0, 150, 100)

#Layer 2
z2		=	depth[1:4]
alphaR	=	cb*np.array([1,1])
betaR	=	0.0*np.array([1,1])
rho		=	pb*np.array([1,1])
alphaI	=	ab*np.array([1,1])	  
betaI	=	0.0*np.array([1,1])
ssp2 = SSPraw(z2, alphaR, betaR, rho, alphaI, betaI)

#	Sound-speed layer specifications
raw = [ssp1, ssp2]
NMedia		=	2
Opt			=	'CVW'	
N			=	[0, 0]	
sigma		=	[0, 0]	
raw[0]
ssp = SSP(raw, depth, 2, Opt, N, sigma)
ssp.make_sspf()


hs = HS(alphaR=cb, betaR=0, rho = pb, alphaI=ab, betaI=0)
Opt = 'A~'
bottom = BotBndry(Opt, hs)
top = TopBndry('CVW')
bdy = Bndry(top, bottom)

plt.title('Sound speed profile')
plt.plot(ssp.sspf(np.linspace(0, 150, 100)), np.linspace(0, 150, 100))
plt.gca().invert_yaxis()
plt.show()

cInt = Empty()
cInt.High = 1e9
cInt.Low = 1400
RMax = max(X)
freq = 50
write_env('py_env.env', 'SCOOTER', 'Pekeris profile', freq, ssp, bdy, pos, [], cInt, RMax)
  
s = Source([sd])
ran =  np.arange(.1,10, 10/1e3)
depth = np.array([]) # doesn't matter for scooter
r = Dom(ran, depth)

pos = Pos(s, r)
options = {'scooter': True}
write_fieldflp('py_env', 'RP', pos, **options)

system("scooter.exe py_env")
system("fields.exe py_env")
[x,x,freq,x,Pos1,pressure]= read_shd('py_env.shd')
pressure = abs(pressure)
pressure = np.squeeze(pressure)
pressure = 10*np.log10(pressure/np.max(pressure))
levs = np.linspace(np.min(pressure), np.max(pressure), 20)
plt.title("greens function contour (relative power DB scale)")
plt.contourf(Pos1.r.range, Pos1.r.depth,pressure[:,:],levels=levs)
plt.xlabel("Range (m)")
plt.ylabel("Depth (m)")
plt.gca().invert_yaxis()
plt.show()

