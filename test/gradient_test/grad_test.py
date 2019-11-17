import numpy as np
import sys
sys.path.append("../")
from os import system
from matplotlib import pyplot as plt
from pyat.pyat.env import *
from pyat.pyat.mode_grad import ModeGrad
from pyat.pyat.readwrite import *
from swellex_helpers.CTD import read_ctd as rc


"""
Use CTD data and a swellex environment to compute the gradient of the mode
with respect to the strength of the strongest EOFs
"""


ctd_dir = '/home/hunter/data/swellex/ctds'
zs, ssps = rc.parse_swellex_ctds(ctd_dir)
u,s,array_dc = rc.get_eofs(zs, ssps)
ur, sr = rc.get_res_matrix(u,s, 1e-2)
num_depths, num_eofs = np.shape(ur)
# extend it to full water column
full_ur = np.zeros((int(216.5*2+1), num_eofs))
full_ur[0:num_depths,:] = ur

sd	=	20
rd = [94.125, 99.755, 105.38, 111.00, 116.62, 122.25, 127.88, 139.12, 144.74, 150.38, 155.99, 161.62, 167.26, 172.88, 178.49, 184.12, 189.76, 195.38, 200.99, 206.62, 212.25]
rr	=	2.5

Z = np.arange(0, 150.5, .5)
X = np.arange(0, 10, .01)


s = Source(sd)
r = Dom(X, Z)
pos = Pos(s,r)
pos.s.depth	= [sd]
pos.r.depth	 = rd
pos.r.range		=	X
pos.Nsd = 1
pos.Nrd = len(rd)
envfil = 's9_env'

[TitleEnv, freq, ssp, bdry, pos1, beam, cint, RMax] = read_env(envfil, 'KRAKEN')
write_env('grad_test.env', 'KRAKEN', TitleEnv, freq, ssp, bdry, pos, beam, cint, RMax)
envfil = 'grad_test'

plt.plot(ssp.raw[0].alphaR)
plt.show()

z1, z2, z3 = ssp.raw[0].z, ssp.raw[1].z, ssp.raw[2].z
freq = 50

E1 = full_ur
E2 = np.identity(len(z2))
E3 = np.identity(len(z3))
dp1 = np.array([.5,.5,.5,.5,.5,.5,.5,.5,.5]).reshape(9,1)
dp2 = np.zeros((len(z2), 1))
dp3 = np.zeros((len(z3), 1))
E = [E1, E2]
plt.plot(E1@dp1)
plt.show()
dp = [dp1, dp2]
mg = ModeGrad(envfil, dp, E)
modegrad, modes0, modes1 = mg.compute_grad()
num_modes = np.shape(modes0.phi)[1]
modegrad = modegrad[0]
dpsi = sum([x*dp[0][i] for x,i in zip(modegrad, range(len(dp[0])))])
plt.plot(modes1.phi[:,0] - modes0.phi[:,0])
plt.plot(dpsi[:,0])
plt.show()
plt.plot(modes0.phi[:,1])
plt.plot(modes1.phi[:,1])
plt.show()

# test inverse thing
inv_maybe = modes0.phi.T@modes0.phi
plt.imshow(inv_maybe.real)
plt.show()
print(inv_maybe)

