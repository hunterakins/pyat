import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d



class Source:
    def __init__(self, depth):
        self.depth = depth

class Dom:
    def __init__(self, ran, depth):
        self.range = ran # in km
        self.depth = depth

class Pos:
    def __init__(self, Source, Dom):
        self.s = Source
        self.r = Dom


class SSPraw:
    def __init__(self,z, alphaR, betaR, rho, alphaI, betaI):
        self.z = z
        self.alphaR = alphaR # sound speed in array
        self.betaR = betaR # not sure
        self.rho = rho 
        self.alphaI = alphaI # atten. in array (alpha (z))
        self.betaI = betaI     
        self.sspf = interp1d(self.z, self.alphaR)

class SSP:
    def __init__(self, raw, depth, NMedia, Opt=None, N=None, sigma=None):
        self.NMedia	= NMedia # number of media layers
        self.Opt = Opt # option: not sure what this is for
        self.N			=	[0, 0]	 # not sure
        self.sigma		=	[0, 0]	 #not sure
        self.depth		=	depth # depth array for layers
        self.raw = raw # list of raw ssp profile
        if NMedia > 1:
            self.ssp_vals = np.zeros((len(self.depth), 1))
            j = 0
            for i in range(NMedia-1):
                layer_depth = self.raw[i].z[-1] # get max depth
                self.sspf = lambda z: np.piecewise(z, [z<layer_depth, z>=layer_depth], [lambda z: self.raw[i].sspf(z), lambda z: self.raw[i+1].sspf(z)])
        else:
            self.ssp_vals = self.raw[0].alphaR
            self.sspf = self.raw[0].sspf

class HS:
    def __init__(self, alphaR=np.array([]), betaR=np.array([]), rho=np.array([]), alphaI=np.array([]), betaI=np.array([])):
        self.alphaR = alphaR
        self.betaR = betaR
        self.rho = rho
        self.alphaI = alphaI
        self.betaI = betaI
        


class BotBndry:
    def __init__(self, Opt, Hs, depth=[], betaI=[0], ):
        self.Opt = Opt # 'A' for analytic or 'CVW' for interpolated ssp
        self.hs = Hs

class TopBndry:
    def __init__(self, Opt):
        self.Opt = Opt

class Bndry:
    def __init__(self, top, bot):
        self.Top = top
        self.Bot = bot


        

class Env:
    def __init__(self, SSP, bndy,dens_array=None):
        self.SSP_array = SSP_array
        self.SSP_depths = SSP.depth
        self.SSP_vals = SSP_array[:,1]
        self.sspf = interp1d(self.SSP_depths, self.SSP)
        self.dens_array = dens_array
        if self.dens_array != None:
            self.rho_depths = self.dens_array[:,0]
            self.rho_vals = self.dens_array[:,1]
            self.rho_f = interp1d(self.rho_depths, self.rho_vals)

class Field:
    def __init__(self, Pos, greens_mat, omega):
        self.Pos = Pos
        self.greens = greens_mat
        self.omega = omega # a list of frequencies


class KernInput:
    def __init__(self, Field_r, Field_s, env):
        self.gr = Field_r.greens_mat
        self.num_rcvrs = len(self.gr)
        self.gs = np.flip(Fields.greens_mat, axis=1) # move source off axis
        self.env = env
        self.Pos = self.Field_r.Pos
        
    



         
