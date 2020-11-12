import numpy as np
from copy import deepcopy
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap as LCM
from scipy.interpolate import interp1d


'''
Objects for interfacing with Michael Porter's models. 
Designed in analog with the types of objects called by read_shd, write_fieldflp, write_env, etc.
'''


class Source:
    def __init__(self, depth):
        self.depth = depth
        self.x = None
        self.y = None

class Dom:
    def __init__(self, ran, depth, offsets=None):
        self.range = ran # in km
        self.depth = depth
        if offsets is not None:
            self.offsets = offsets

class Pos:
    def __init__(self, Source, Dom):
        self.s = Source
        self.r = Dom

class Ice:
    def __init__(self, BumDen, eta, xi):
        self.BumDen=  BumDen
        self.eta = eta
        self.xi = xi

class SSPraw:
    def __init__(self,z, alphaR, betaR, rho, alphaI, betaI):
        """
        Input:
        z - numpy 1d array
            Depths the ssp is taken at
        alphaR - numpy 1d array
            ssp vals
        betaR - numpy 1d array
            shear speeds
        rho - numpy 1d array
            density vals at the points
        alphaI - numpy 1d array
            attenuation (p-wave)
        betaI - numpy 1d array
            shear attenuation
        """
        self.z = z
        self.alphaR = alphaR # sound speed in array
        self.betaR = betaR # shear wave speed
        self.rho = rho 
        self.alphaI = alphaI # atten. in array (alpha (z))
        self.betaI = betaI # shear attenuation

    def make_sspf(self):
        self.sspf = interp1d(self.z, self.alphaR)

    def interp_all(self):
        self.betaI_f = interp1d(self.z, self.betaI)
        self.betaR_f = interp1d(self.z, self.betaR)
        self.rho_f = interp1d(self.z, self.rho)
        self.alphaR_f = interp1d(self.z, self.alphaR)
        self.alphaI_f = interp1d(self.z, self.alphaI)
        return self.alphaR_f, self.betaR_f, self.rho_f, self.alphaI_f, self.betaI_f

class SSP:
    def __init__(self, raw, depth, NMedia, Opt=None, N=None, sigma=0, ranges=0):
        self.NMedia	= NMedia # number of media layers
        self.Opt = Opt # option: not sure what this is for
        self.N			=	N	 # array with num points in each layer
        self.sigma		=	sigma	 # not sure
        self.depth		=	depth # depth array for layers
        self.ranges = ranges
        self.raw = raw # list of raw ssp profile
        self.sspf = None

    def make_sspf(self):
        for raw in self.raw:
            raw.make_sspf()
        if self.NMedia > 1:
            self.ssp_vals = np.zeros((len(self.depth), 1))
            if self.NMedia == 2:
                layer_depth = self.raw[0].z[-1]
                self.sspf = lambda z: np.piecewise(z, [z<layer_depth, z>=layer_depth], [lambda z: self.raw[0].sspf(z), lambda z: self.raw[1].sspf(z)])
            elif self.NMedia == 3:
                layer_depth_one = self.raw[0].z[-1]
                layer_depth_two = self.raw[1].z[-1]
                func1 = lambda z: np.piecewise(z, [z<layer_depth_one, z>= layer_depth_one], [lambda z: self.raw[0].sspf(z), lambda z: self.raw[1].sspf(z)])
                self.sspf = lambda z: np.piecewise(z, [z<layer_depth_two, z>= layer_depth_two], [lambda z: func1(z), lambda z: self.raw[2].sspf(z)])
            else:
                raise ValueError("too many layers")
        else:
            self.ssp_vals = self.raw[0].alphaR
            self.sspf = self.raw[0].sspf

class HS:
    def __init__(self, alphaR=np.array([]), betaR=np.array([]), rho=np.array([]), alphaI=np.array([]), betaI=np.array([])):
        self.alphaR = np.array(alphaR)
        self.betaR = np.array(betaR)
        self.rho = np.array(rho)
        self.alphaI = np.array(alphaI)
        self.betaI = np.array(betaI)
        
class BotBndry:
    def __init__(self, Opt, Hs, depth=[], betaI=[0], ):
        self.Opt = Opt # 'A' for analytic or 'CVW' for interpolated ssp
        self.hs = Hs

class TopBndry:
    def __init__(self, Opt, depth=[]):
        self.Opt = Opt
        self.cp = None
        self.cs = None
        self.rho = None 

class Bndry:
    def __init__(self, top, bot):
        self.Top = top
        self.Bot = bot

class Box:
    def __init__(self, z, r):
        self.r = r
        self.z = z

class Beam:
    def __init__(self, RunType=None, Type=None,Nbeams=None, Ibeam=None, Nrays=None, alpha=None, deltas=None, box=None, epmult=None, rloop=None, Ibwin=None, Nimage = None):
        self.RunType = RunType
        self.Type = Type
        self.Nbeams =  Nbeams
        self.Ibeam  =  Ibeam 
        self.Nrays  =  Nrays
        self.alpha  =  alpha
        self.deltas =  deltas
        self.box    =  box
        self.epmult =  epmult
        self.rloop  =  rloop
        self.Ibwin = Ibwin
        self.Nimage = None

class cInt:
    def __init__(self, low, high):
        self.Low = low
        self.High = high

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
    def __init__(self, **kwargs):
        self.pos = kwargs['pos']
        self.greens = kwargs['pressure']
        self.freq = kwargs['freq'] # a list of frequencies

    def contour(self):
        shape = np.shape(self.greens)
        if len(shape) > 2:
            raise ValueError("Field object is too high dimensional, try multi_contour")
        pos = self.pos 
        ax = plt.contourf(pos.r.range, -pos.r.depth, np.log10(abs(self.greens.real)))
        return ax

    def multi_contour(self):
        print('not implemented, maybe one day...')

class Modes:
    def __init__(self, **kwargs):
        self.M = kwargs['M']
        self.k = kwargs['modes_k']
        self.z = np.array(kwargs['z'])
        self.phi = kwargs['modes_phi']
        self.top = kwargs['top']
        self.bot = kwargs['bot']
        self.N = kwargs['N'] 
        self.Nfreq = kwargs['Nfreq']
        self.Nmedia = kwargs['Nmedia']
        self.depth = kwargs['depth']
        self.rho = kwargs['rho']
        self.freqvec = kwargs['freqVec']
        self.init_dict = kwargs
        self.num_modes = self.M # easier to remember

    def get_excited_modes(self, sd, threshold):
        '''
        return an array of modes that are excited by a source at sd meters depth
        threshold 
        
        also populate some structures in moded caled excited_phi and excited_k

        '''
        if sd not in self.z:
            raise ValueError("sd not in the depth array, are you sure it's the right depth you're passing?")
        depth_ind = [i for i in range(len(self.z)) if self.z[i] == sd][0]
        vals = self.phi[depth_ind,:]
        const = np.max(abs(vals)) 
        filtered_inds = [i for i in range(len(self.k)) if abs(self.phi[depth_ind, i]) / const > threshold]
        self.excited_phi = self.phi[:, filtered_inds]
        self.excited_k = self.k[filtered_inds]
        return self.excited_phi, self.excited_k

    def get_source_depth_ind(self, sd):
        if sd not in self.z:
            raise ValueError("sd not in the depth array, are you sure it's the right depth you're passing?")
        else:
            sind = [i for i in range(len(self.z)) if self.z[i] == sd][0]
            self.sind = sind
            return  sind

    def remove_source_pos(self, sd):
        """
        Take the source at sd from the mode matrix
        Initiate a new attribute to hold the source
        modal value
        """
        sind = self.get_source_depth_ind(sd)
        new_pos_len = len(self.z) - 1
        new_phi = np.zeros((new_pos_len, self.num_modes), dtype=self.phi.dtype)
        new_phi[:sind, :] = self.phi[:sind, :]
        new_phi[sind:,:] = self.phi[sind+1:,:]
        self.phi = new_phi
        new_z = np.zeros((new_pos_len), dtype=self.z.dtype)
        new_z[:sind] = self.z[:sind]
        new_z[sind:] = self.z[sind+1:]
        self.z = new_z
        self.source_strength = self.phi[sind,:]
        return 
            
    def get_source_strength(self, sd):
        """
        Get the value of each mode at the source depth sd (meters)
        Initialize new attribute for the source strength
        """
        sind = self.get_source_depth_ind(sd)
        vals = self.phi[sind,:]
        self.source_strength = vals
        return  vals

    def get_receiver_modes(self, zr):
        r_inds = [i for i in range(len(self.z)) if self.z[i] in zr]
        receiver_modes = self.phi[r_inds, :]
        self.receiver_modes = receiver_modes
        return receiver_modes
        
        

    def plot(self):
        figs = []
        if self.M > 5:
            for i in range(self.M):
                fig = plt.figure(i)
                plt.plot(self.phi[:,i], -self.z)
                figs.append(fig)
        return figs

    def __repr__(self):
        return 'Modes object with ' + str(self.M) + ' distinct modes'
      
class Arrival:
    def __init__(self, info_list):
        self.amp = info_list[0]
        self.delay = info_list[1] 
        self.src_ang = info_list[2] 
        self.rec_ang = info_list[3] 
        self.num_top_bnc = info_list[4] 
        self.num_bot_bnc = info_list[5] 
        
class Arrivals:
    def __init__(self, arrival_list):
        self.arrivals = arrival_list
    
    def plot_cir(self, vals=None):
        """
        Plot the channel impulse response
        Vals provides the option to pass in a set of arrivals to use fo calibrating the y axis
        If vals=None, I use the arrivals in the object to calibrate the y-axis
        I don't create a figure
        I return ymin, ymax, vals so that the values from this arrival set can be used in a more global
        scope to set the y-axis limits (useful for creating videos comprised of snapshots of a bunch
        of CIR's). For example, I want to plot a bunch of Arrivals on the same plot, or combine a bunch of
        CIRs into a video. Then I need to have the same axis limits for all the plots. 
        I can use default vals (None) for the first plot, then use the ymin, ymax, vals return values
        in a larger scope to standardize the axes. Then for the second, third, ... plots, I pass in vals
        to make sure they all have the same time axis.
        Output - 
            ymin - float
            ymax -float 
                plus/minus 1.5 times max amplitude value
            vals - list
                vals[0] is the times
                vals[1] is the amplitude vals
        """
        arrival_list = self.arrivals
        amps = np.array([x.amp.real for x in arrival_list])
        phase = np.array([np.angle(x.amp) for x in arrival_list])
        times = np.array([x.delay for x in arrival_list])
        """
        Create the time axis
        """
        if type(vals) != type(None):
            t = np.linspace(np.min(vals[0])-2, np.max(vals[0]) +2, 100)
        else:
            t = np.linspace(np.min(times) - 2, np.max(times) + 2, 100)
        zeros = np.zeros(t.size)
        plt.scatter(t, zeros, s=6)
        plt.stem(times, amps.real, markerfmt=' ', basefmt=' ',use_line_collection=True)
        plt.scatter(times, amps,s=36,c='r')
        scale = 1.5*np.max(abs(amps)) # set yaxis scale
        ymin = - scale
        ymax = scale
        vals = [times, amps]
        return ymin, ymax, vals

    def plot_ellipse(self, vals=None):
        """
        Produce the travel time ellipse with colorbar for amplitude
        Input - 
        self - contains all the arrival info
        vals - optional
            allows me to use a reference arrival set to calibrate the plot axes
            vals[0] = times
            vals[0] = amps
        """
        arrival_list = self.arrivals
        amps = np.array([x.amp.real for x in arrival_list])
        times = np.array([x.delay for x in arrival_list])
        angles = np.array([x.rec_ang for x in arrival_list])
        """
        Create the time axis
        """
        if type(vals) != type(None):
            tmin, tmax = np.min(vals[0]), np.max(vals[0])
        else:
            tmin, tmax = np.min(times), np.max(times)
        amps = amps / np.max(abs(amps))
        cvals = np.log10(abs(amps))
        """ Make a nice cmap """ 
        oran = cm.get_cmap('hsv')
        oran = oran(np.linspace(.13, 0, 10))
        mymap = LCM(oran, 'mymap')
        """ Configure the x axis """
        plt.scatter([tmin, tmax], [0,0], s=0)
        plt.scatter(times, angles, s=25, c=cvals, cmap=mymap)
        plt.colorbar()
        vals = [times, amps]
        ymin, ymax = -1.5*np.max(abs(angles)), 1.5*np.max(abs(angles))
        return ymin, ymax, vals

class KernInput:
    def __init__(self, Field_r, Field_s, env):
        self.gr = Field_r.greens_mat
        self.num_rcvrs = len(self.gr)
        self.gs = np.flip(Fields.greens_mat, axis=1) # move source off axis
        self.env = env
        self.Pos = self.Field_r.Pos
        
    



         
