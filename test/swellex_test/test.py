import numpy as np
import sys
sys.path.append("../")
from os import system
from matplotlib import pyplot as plt
from pyat.pyat.env import *
from pyat.pyat.readwrite import *





"""
Environment specification here
"""

pw		=	1.0
aw		=	0.0
Z1 = 216.5
Z2 = 240.0
depth = [0.0, Z1, Z2]
# Layer 1
z1 = np.array([0.0, 0.83, 1.39, 2.13, 2.95, 3.65, 4.53, 5.44, 6.15, 6.8, 7.73, 8.69, 9.45, 10.12, 11.0, 12.0, 12.79, 13.53, 14.29, 15.27, 16.1, 16.85, 17.87, 19.03, 19.78, 20.33, 21.16, 22.17, 23.25, 24.5, 25.77, 26.97, 28.12, 29.1, 29.5, 29.73, 29.9, 30.27, 30.59, 30.98, 31.19, 31.31, 31.44, 31.81, 31.94, 32.08, 32.33, 32.55, 32.71, 32.9, 33.41, 34.17, 35.03, 35.89, 36.78, 37.82, 38.9, 39.81, 40.48, 41.16, 42.24, 43.42, 44.51, 45.65, 46.72, 47.72, 48.72, 49.65, 50.6, 51.08, 51.48, 51.85, 51.98, 52.68, 53.81, 54.86, 55.96, 57.07, 58.15, 59.0, 59.73, 60.6, 61.77, 62.85, 63.85, 64.82, 66.0, 66.96, 67.83, 68.74, 69.74, 70.71, 71.69, 72.65, 73.7, 74.72, 75.93, 77.04, 78.0, 78.83, 79.65, 80.5, 81.32, 81.94, 82.69, 83.79, 84.95, 86.11, 87.25, 88.26, 89.06, 89.88, 90.7, 91.39, 92.02, 92.73, 93.35, 93.85, 94.59, 95.56, 96.08, 96.4, 96.69, 96.88, 97.34, 97.79, 98.42, 98.95, 99.06, 99.28, 99.59, 100.5, 101.39, 102.16, 102.88, 103.66, 104.46, 104.97, 105.47, 106.24, 107.04, 107.6, 108.18, 108.84, 109.53, 109.98, 110.23, 110.71, 111.21, 111.44, 111.68, 112.09, 112.38, 112.75, 113.2, 113.61, 113.93, 114.22, 114.8, 115.94, 117.19, 118.56, 119.98, 121.38, 122.57, 123.74, 125.02, 126.39, 127.67, 129.04, 130.4, 131.52, 132.62, 133.9, 135.11, 136.24, 137.41, 138.69, 139.93, 140.97, 142.04, 143.16, 144.22, 145.22, 146.23, 147.28, 148.25, 149.29, 150.45, 151.46, 152.52, 153.76, 154.85, 155.8, 156.82, 157.94, 159.05, 160.21, 161.33, 162.35, 163.46, 164.52, 164.75, 168.17, 169.2, 170.18, 171.13, 172.08, 172.96, 173.84, 174.68, 175.51, 176.39, 177.13, 177.9, 178.75, 179.73, 180.38, 181.09, 181.78, 182.41, 183.02, 183.57, 184.15, 184.9, 185.59, 186.25, 187.01, 187.95, 188.96, 189.98, 191.04, 192.18, 193.13, 193.9, 194.43, 195.05, 195.71, 196.17, 196.67, 197.27, 197.92, 198.47, 198.97, 199.32, 199.59, 199.95, 200.16, 208.35, 216.5])
alphaR	=	np.array([1521.94, 1521.95, 1521.94, 1521.9, 1521.84, 1521.74, 1521.53, 1521.18, 1520.84, 1520.51, 1520.04, 1519.58, 1519.13, 1518.59, 1517.45, 1515.25, 1512.28, 1509.48, 1507.26, 1505.05, 1503.4, 1502.35, 1501.27, 1500.28, 1499.89, 1499.64, 1499.22, 1498.59, 1497.96, 1497.13, 1496.52, 1496.12, 1495.77, 1495.49, 1495.4, 1495.34, 1495.3, 1495.21, 1495.12, 1495.02, 1494.95, 1494.91, 1494.87, 1494.7, 1494.63, 1494.56, 1494.4, 1494.27, 1494.21, 1494.14, 1493.93, 1493.67, 1493.44, 1493.25, 1493.1, 1492.92, 1492.68, 1492.54, 1492.45, 1492.37, 1492.27, 1492.18, 1492.13, 1492.08, 1492.03, 1491.98, 1491.93, 1491.89, 1491.84, 1491.82, 1491.79, 1491.77, 1491.76, 1491.69, 1491.54, 1491.37, 1491.23, 1491.13, 1491.04, 1490.97, 1490.91, 1490.84, 1490.74, 1490.66, 1490.59, 1490.51, 1490.42, 1490.35, 1490.3, 1490.24, 1490.17, 1490.07, 1489.93, 1489.85, 1489.79, 1489.73, 1489.68, 1489.66, 1489.63, 1489.61, 1489.61, 1489.62, 1489.62, 1489.63, 1489.63, 1489.63, 1489.62, 1489.59, 1489.55, 1489.51, 1489.46, 1489.41, 1489.33, 1489.24, 1489.16, 1489.14, 1489.13, 1489.12, 1489.07, 1488.96, 1488.89, 1488.83, 1488.78, 1488.75, 1488.65, 1488.59, 1488.55, 1488.51, 1488.5, 1488.49, 1488.46, 1488.33, 1488.25, 1488.2, 1488.18, 1488.17, 1488.16, 1488.15, 1488.13, 1488.08, 1488.05, 1488.03, 1488.0, 1487.99, 1487.97, 1487.96, 1487.95, 1487.92, 1487.87, 1487.84, 1487.81, 1487.77, 1487.75, 1487.73, 1487.71, 1487.7, 1487.69, 1487.69, 1487.7, 1487.73, 1487.78, 1487.83, 1487.91, 1488.02, 1488.09, 1488.14, 1488.19, 1488.28, 1488.4, 1488.54, 1488.67, 1488.74, 1488.81, 1488.87, 1488.92, 1488.95, 1488.97, 1488.96, 1488.95, 1488.93, 1488.92, 1488.91, 1488.91, 1488.91, 1488.91, 1488.89, 1488.89, 1488.89, 1488.9, 1488.9, 1488.92, 1488.93, 1488.95, 1488.96, 1488.98, 1488.99, 1489.0, 1489.0, 1488.99, 1488.96, 1488.92, 1488.9, 1488.88, 1488.53, 1488.47, 1488.42, 1488.44, 1488.43, 1488.45, 1488.46, 1488.48, 1488.48, 1488.5, 1488.51, 1488.53, 1488.57, 1488.55, 1488.6, 1488.61, 1488.61, 1488.55, 1488.55, 1488.54, 1488.54, 1488.55, 1488.53, 1488.5, 1488.51, 1488.53, 1488.51, 1488.47, 1488.4, 1488.35, 1488.34, 1488.32, 1488.31, 1488.27, 1488.25, 1488.25, 1488.26, 1488.27, 1488.29, 1488.3, 1488.3, 1488.32, 1488.31, 1488.31, 1488.26, 1488.26, 1488.26])
betaR	=	0.0*np.ones(z1.shape) # this is shear speed
rho		=	pw*np.ones(z1.shape) # density
alphaI	=	aw*np.ones(z1.shape) # attneuation
betaI	=	0.0*np.ones(z1.shape) # shear attenuation
ssp1 = SSPraw(z1, alphaR, betaR, rho, alphaI, betaI)

# layer 2
# sediment layer parameters

z2 = np.array([Z1, Z2])
alphar2 = np.array([1572.37, 1593.02])
betaR2 = np.zeros(z2.shape)
rho2 = np.array([1.76, 1.76])
alphaI2 = np.array([0.2, 0.2])
ssp2 = SSPraw(z2, alphar2, betaR2, rho2, alphaI2, betaI)

# clay layer is implemented as a halfspace here (traditionally it is included as an 800 m thick layer with basalt halfspace but there is no propagation sensitivity below ~ 240 m)
cb = 1800.0
pb = 2.0
ab = 0.06


raw = [ssp1, ssp2]
NMedia		=	2 # number of layers
Opt			=	'CVW'	 # C means linearly interpolate sound speed V means vacuum surface (pressure release) W means  that attenuation is specified in db per lambda (see Kraken documentation on env files)
N			=	[z1.size, z2.size]
sigma		=	None
ssp = SSP(raw, depth, NMedia, Opt, N, sigma)


# Acoustic halfpsace
hs = HS(alphaR=cb, betaR=0, rho = pb, alphaI=ab, betaI=0)
Opt = 'A~' #this means acoustic halspace below
bottom = BotBndry(Opt, hs)
top = TopBndry('CVW')
bdy = Bndry(top, bottom)


class Empty:
    def __init__(self):
        return

# this determines the lower and upper modal phase speed cutoffs
# read Chapter 2. of Computational Ocean Acoustics for more information
cInt = Empty()
cInt.High = cb
cInt.Low = 0 # compute automatically


"""
Geometry is here
"""
freq = 200


sd=	50.0
rd = [94.125, 99.755, 105.38, 111.00, 116.62, 122.25, 127.88, 139.12, 144.74, 150.38, 155.99, 161.62, 167.26, 172.88, 178.49, 184.12, 189.76, 195.38, 200.99, 206.62, 212.25]
rr	=	2.5
s = Source(sd)

Z = np.arange(0, Z1, .5)
X = np.arange(0, 10, .01) # in km
RMax = max(X)
r = Dom(X, Z)
pos = Pos(s,r)
pos.s.depth	= [sd]
pos.r.depth	 = Z
pos.r.range		=	X
pos.Nsd = 1
pos.Nrd = len(rd)

write_env('swellex_env.env', 'KRAKEN', 'SWellEx environment', freq, ssp, bdy, pos, [], cInt, RMax)
  

s = Source([sd])
ran =  np.arange(10/1e3, 10, 10/1e3)
depth = np.arange(0,Z1,1)
r = Dom(ran, depth)

pos = Pos(s, r)

write_fieldflp('swellex_env', 'R', pos)
system("krakenc.exe swellex_env")
fname = 'swellex_env.mod'
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
system("field.exe swellex_env")
[x,x,x,x,Pos1,pressure]= read_shd('swellex_env.shd')
#pressure = np.squeeze(pressure)
pressure = abs(pressure)
pressure = 20*np.log10(pressure / np.max(pressure))
print(np.shape(pressure))
levs = np.linspace(np.max(pressure)-60, np.max(pressure), 20)
plt.contourf(Pos1.r.range, Pos1.r.depth,(pressure[0, 0,:,:]),levels=levs)
plt.gca().invert_yaxis()
plt.xlabel('Range (km)')
plt.ylabel('Depth (m)')
plt.show()

