import numpy as np
import sys
sys.path.append("../")
from os import system
from matplotlib import pyplot as plt
from pyat.pyat.env import *
from pyat.pyat.mode_grad import ModeGrad
from pyat.pyat.readwrite import *



sd	=	20
rd = [94.125, 99.755, 105.38, 111.00, 116.62, 122.25, 127.88, 139.12, 144.74, 150.38, 155.99, 161.62, 167.26, 172.88, 178.49, 184.12, 189.76, 195.38, 200.99, 206.62, 212.25]
rr	=	2.5

Z = np.arange(0, 150.5, .5)
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
print(alphaR)
betaR	=	0.0*np.array([1,1])		
rho		=	pw*np.array([1,1])		
alphaI	=	aw*np.array([1,1])		
betaI	=	0.0*np.array([1,1])

ssp1 = SSPraw(z1, alphaR, betaR, rho, alphaI, betaI)
print(ssp1.alphaR)
arr = np.linspace(0, 150, 100)

#Layer 2
z2		=	depth[1:]
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

hs = HS()
Opt = 'A~'
bottom = BotBndry(Opt, hs)
top = TopBndry('CVW')
bdy = Bndry(top, bottom)

plt.plot(ssp.sspf(np.linspace(0, 150, 100)))
plt.show()

class Empty:
    def __init__(self):
        return

cInt = Empty()
cInt.High = 1e9
cInt.Low = 1400
RMax = max(X)
freq = 50
E1 = np.array([[1,2], [1, -2]])
print(E1)
E2 = np.identity(2)
dp1 = np.array([1,-1]).reshape(2,1)
dp2 = np.array([-10,-10]).reshape(2,1)
E = [E1, E2]
dp = [dp1, dp2]
mg = ModeGrad('py_env', dp, E)
modegrad, modes0, modes1 = mg.compute_grad()
plt.plot(modegrad[:,0])
plt.show()
plt.plot(modes0.phi[:,0])
plt.show()
plt.plot(modes1.phi[:,0])
plt.show()

write_env('py_env.env', 'KRAKEN', 'Pekeris profile', freq, ssp, bdy, pos, [], cInt, RMax)
  

s = Source([sd])
ran =  np.arange(0,10, 10/1e3)
depth = np.arange(0,150,1)
r = Dom(ran, depth)
pos = Pos(s, r)


