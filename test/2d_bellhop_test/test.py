import numpy as np
import sys
sys.path.append("../")
from os import system
from matplotlib import pyplot as plt
from pyat.pyat.env import *
from pyat.pyat.readwrite import *

"""
Run a bellhop test on a simple range-dependent environment Munk profile
Goal is look at contour map of green's function
"""


"""
Basic source receiver parameters
"""

sd	=	1000 # meters source depth
rr	=	2.5

"""
Grid for simulation (each point in the grid is a receiver)
"""
Z = np.linspace(1, 5000, 201) # in meters
X = np.arange(1, 100, .1) # in kilometers



"""
Initialize pyat objects with the parameters configured above
"""
s = Source(sd)
r = Dom(X, Z)
pos = Pos(s,r)
pos.s.depth	= [sd] # something weird about how this populates so have to manually do it
pos.r.depth	 = Z
pos.r.range		=	X
pos.Nsd = 1 # number of source ranges
pos.Nrd = len(X) # number of receiver ranges

"""
Set up environment
"""
cw		=	1500 # default sound speed
pw		=	1 # density (kg/m^3
aw		=	0 # atten in water

depth = [0, 5000] # layer depths (surface and bottom)
z1 = [0.0,  200.0,  250.0,  400.0,  600.0,  800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0] # SSP profile depths

alphaR = [1548.52,1530.29,1526.69,1517.78,1509.49,1504.30,1501.38,1500.14,1500.12,1501.02,1502.57,1504.62,1507.02,1509.69,1512.55,1515.56,1518.67,1521.85,1525.10,1528.38,1531.70,1535.04,1538.39,1541.76,1545.14,1548.52,1551.91] # SSP profile values
betaR	=	0.0*np.ones((len(z1))) # shear speed at each depth
rho		=	pw*np.ones((len(z1))) # rho at each depth
alphaI	=	aw*np.ones((len(z1))) # atten at each depth
betaI	=	0.0*np.ones((len(z1))) # shear atten at each depth

# load the SSP into a 2D numpy array for write_env
raw = np.zeros((len(alphaR), 2))
raw[:,0] = alphaR
"""
Create a second SSP by simply adding 10 to every value in the original Munk SSP
"""
raw[:,1] = [x + 10 for x in alphaR] 
ssp_ranges = np.array([0,110]) # in km

"""
Write the environment to a ssp file (required only for range dependent bellhop
"""
write_ssp('py_env', raw, ssp_ranges)

# Create SSPraw object out of the first column
ssp1 = SSPraw(z1, alphaR, betaR, rho, alphaI, betaI)
raw = [ssp1] # make it a list type
NMedia		=	1 
Opt			=	'QVW'	 # see AT docs for details on env options
N			=	[0, 0]	
sigma		=	[0, 0]	
# SSP object contains all the info for the first layer (the water layer)
ssp = SSP(raw, depth, NMedia, Opt, N, sigma, 2)

"""
Set up second layer (the sediment bottom
"""
alphaR = 1600 # p wave speed in sediment
betaR = 0 # no shear wave 
alphaI = .5 # p wave atten
betaI =0 # s wave atten
rhob = 1.6 # density

"""
Halfspace object treats the second layer as an infinite halfspace
with bottom sound speed parameters
"""
hs = HS(alphaR, betaR, rhob, alphaI, betaI)
Opt = 'A' # analytical solution in the bottom matches boundary conditions in water column
bottom = BotBndry(Opt, hs)
top = TopBndry('QVW')
bdy = Bndry(top, bottom) # package top and bottom together


# min and max ssp (only used for modal calculations)
low = 1400
high = 1e9
cInt = cInt(low, high)
RMax = max(X)
freq = 3000

# Beam params
run_type = 'I' # incoherently sum beams, see AT docs for more info
nbeams = 100
alpha = np.linspace(-20,20, nbeams) # min and max launch angle -20 degrees to 20 degrees
box = Box(5500, 100) #bound the region you let the beams go, depth in meters and range in km
deltas=0 # length step of ray trace, 0 means automatically choose
beam = Beam(RunType=run_type, Nbeams=nbeams, alpha=alpha,box=box,deltas=deltas) # package

"""
Write the env file for BELLHOP
"""
write_env('py_env.env', 'BELLHOP', 'Pekeris profile', freq, ssp, bdy, pos, beam, cInt, RMax)

"""
Run bellhop
"""
system("bellhop.exe py_env")

"""
Read in shade file to numpy array
"""
[x,x,x,x,ppos, p] = read_shd("py_env.shd")
p = np.squeeze(p)
p = abs(p)
p = 10*np.log10(p/np.max(p)) # convert to relative power
levs = np.linspace(-30, 0, 20) # set the levels of the contour map

"""
Plot it
"""
plt.title('Contour for 2-D Munk profile')
plt.contourf(np.squeeze(p), levels=levs) 
plt.gca().invert_yaxis()
plt.show()


