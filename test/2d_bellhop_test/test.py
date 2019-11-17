import numpy as np
import sys
sys.path.append("../")
from os import system
from matplotlib import pyplot as plt
from pyat.pyat.env import *
from pyat.pyat.readwrite import *



sd	=	1000
#rd = [94.125, 99.755, 105.38, 111.00, 116.62, 122.25, 127.88, 139.12, 144.74, 150.38, 155.99, 161.62, 167.26, 172.88, 178.49, 184.12, 189.76, 195.38, 200.99, 206.62, 212.25]
rr	=	2.5

Z = np.linspace(1, 5000, 201)
X = np.arange(1, 100, .1)

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


depth = [0, 5000]
z1 = [0.0,  200.0,  250.0,  400.0,  600.0,  800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0]

alphaR = [1548.52,1530.29,1526.69,1517.78,1509.49,1504.30,1501.38,1500.14,1500.12,1501.02,1502.57,1504.62,1507.02,1509.69,1512.55,1515.56,1518.67,1521.85,1525.10,1528.38,1531.70,1535.04,1538.39,1541.76,1545.14,1548.52,1551.91]
betaR	=	0.0*np.ones((len(z1)))
rho		=	pw*np.ones((len(z1)))
alphaI	=	aw*np.ones((len(z1)))
betaI	=	0.0*np.ones((len(z1)))


raw = np.zeros((len(alphaR), 2))
raw[:,0] = alphaR
raw[:,1] = [x + 10 for x in alphaR]
write_ssp('py_env', raw, np.array([0, 110]))



ssp1 = SSPraw(z1, alphaR, betaR, rho, alphaI, betaI)


#	Sound-speed layer specifications
raw = [ssp1]
NMedia		=	1
Opt			=	'QVW'	
N			=	[0, 0]	
sigma		=	[0, 0]	
raw[0]
ssp = SSP(raw, depth, NMedia, Opt, N, sigma, 2)


# Layer 2
alphaR = 1600 # p wave speed in sediment
betaR = 0 # no shear wave
alphaI = .5 # p wave atten
betaI =0 # s wave atten
rhob = 1600

hs = HS(alphaR, betaR, rhob, alphaI, betaI)
Opt = 'A'
bottom = BotBndry(Opt, hs)
top = TopBndry('QVW')
bdy = Bndry(top, bottom)

#plt.plot(ssp.sspf(np.linspace(0, 4000, 100)))
#plt.show()

class Empty:
    def __init__(self):
        return

low = 1400
high = 1e9
cInt = cInt(low, high)
RMax = max(X)
freq = 3000

# Beam params
run_type = 'I'
nbeams = 100
alpha = np.linspace(-20,20, 100)
box = Box(5500, 100)
deltas=0
beam = Beam(RunType=run_type, Nbeams=nbeams, alpha=alpha,box=box,deltas=deltas)

write_env('py_env.env', 'BELLHOP', 'Pekeris profile', freq, ssp, bdy, pos, beam, cInt, RMax)
  

#s = Source([sd])
#ran =  np.arange(0,10, 10/1e3)
#depth = np.arange(0,150,1)
#r = Dom(ran, depth)


system("bellhop.exe py_env")

[x,x,x,x,ppos, p] = read_shd("py_env.shd")
print(p.shape)
p = np.squeeze(p)
p = abs(p)
p = p[:,:]
p = 10*np.log10(p/np.max(p))
print(p)
levs = np.linspace(-30, 0, 20)
plt.contourf(np.squeeze(p), levels=levs)
plt.gca().invert_yaxis()
plt.show()

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
