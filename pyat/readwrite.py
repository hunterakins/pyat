import numpy as np
from scipy.io import loadmat
import os
from struct import unpack
from .env import Source, Dom, Pos, cInt, Ice, SSPraw, SSP, HS, BotBndry, TopBndry, Bndry, Box, Beam, Modes
import re
from matplotlib import pyplot as plt


""" 
    translation of basic read functions from Mike Porter's 
    Matlab ReadWrite folder

    Currently lacking bellhop support (just need to translate read_bell)
    and bellhop3d (need to translate readsxrx)
"""


def equally_spaced( x ):
    # routine to see if the vector x is composed of roughly equally spaced
    # values
    n  = len( x )
    xtemp = np.linspace( x[ 0 ], x[ -1 ], n )
    delta = abs( x[:] - xtemp[:] )   # the colon converts both to column vectors

    if ( np.max( delta ) < 1e-9 ):
        return  1
    else:
        return 0

# Write a field-parameters file
def write_fieldflp(flpfil, Option, Pos, **kwargs):
    if flpfil[-4:] != '.flp':
        flpfil = flpfil + '.flp'

    if 'scooter' in kwargs.keys():
        with open(flpfil, 'w') as f:
            f.write('\'' + '{:4}'.format(Option) + ' \'' + ' ! Option \r\n')
            f.write('{:5d}'.format(len(Pos.r.range)) + ' \t \t \t \t ! NRR')
            if ( len( Pos.r.range ) > 2) and equally_spaced( Pos.r.range ):
                f.write('\r\n    ' + '{:6f}'.format(Pos.r.range[0]))
                f.write('\r\n    ' + '{:6f}'.format(Pos.r.range[-1]))
            elif len(Pos.r.range) == 1:
                f.write('\r\n    ' + '{:6f}'.format(Pos.r.range[0]) + '  ')
            else:
                for i in range(len(Pos.r.range)):
                    f.write('\r\n    ' + '{:6f}'.format(Pos.r.range[i]) + '  ')
                

            f.write( '/ \t ! RR(1)  ... (km) \r\n' )
        return
    with open(flpfil, 'w') as f:
        f.write('/ ! Title \r\n' )
        f.write('\'' + '{:4}'.format(Option) + ' \'' + ' ! Option \r\n')
        f.write('999999   ! Mlimit (number of modes to include) \r\n' )
        f.write('1   ! NProf \r\n' )
        f.write('0.0 /  ! rProf (km) \r\n' )

        # receiver ranges
        f.write('{:5d}'.format(len(Pos.r.range)) + ' \t \t \t \t ! NRR')

        if ( len( Pos.r.range ) > 2) and equally_spaced( Pos.r.range ):
            f.write('\r\n    ' + '{:6f}'.format(Pos.r.range[0]))
            f.write('\r\n    ' + '{:6f}'.format(Pos.r.range[-1]))
        elif len(Pos.r.range) == 1:
            f.write('\r\n    ' + '{:6f}'.format(Pos.r.range[0]) + '  ')
        else:
            for i in range(len(Pos.r.range)):
                f.write('\r\n    ' + '{:6f}'.format(Pos.r.range[i]) + '  ')
            

        f.write( '/ \t ! RR(1)  ... (km) \r\n' )

        # source depths

        f.write('{:5d}'.format(len(Pos.s.depth))+'\t \t \t \t ! NSD')

        if ( len( Pos.s.depth ) > 2 and equally_spaced( Pos.s.depth ) ):
            f.write('\r\n    '+'{:6f}'.format(Pos.s.depth[ 0 ])+'  ')
            f.write('\r\n    '+'{:6f}'.format(Pos.s.depth[ -1 ])+'  ')
        elif len(Pos.s.depth) == 1:
            f.write( '\r\n    '+ '{:6f}'.format(Pos.s.depth[0])+'  ')
        else:
            for i in range(len(Pos.s.depth)):
                f.write( '\r\n    '+ '{:6f}'.format(Pos.s.depth[i])+'  ')

        f.write( '/ \t ! SD(1)  ... (m) \r\n' )

        # receiver depths

        f.write('{:5}'.format(len(Pos.r.depth)) + ' \t \t \t \t ! NRD')
    
        if ( len( Pos.r.depth ) > 2 and equally_spaced( Pos.r.depth ) ):
            f.write( '\r\n    '+'{:6f}'.format(Pos.r.depth[0])+'  ')
            f.write( '\r\n    '+'{:6f}'.format(Pos.r.depth[-1])+'  ')
        elif (len(Pos.r.depth) == 1):
            f.write( '\r\n    '+'{:6f}'.format(Pos.r.depth[0])+'  ')
        else:
            for i in range(len(Pos.r.depth)):
                f.write( '\r\n    '+'{:6f}'.format(Pos.r.depth[i]))

        f.write( '/ \t ! RD(1)  ... (m) \r\n' )

        # receiver range offsets
        f.write( '{:5}'.format(len( Pos.r.depth ))+' \t \t ! NRR',  )
        if not hasattr(Pos.r, 'offsets'):
            f.write( '\r\n    '+'{:6.2f}'.format(0.0) +'  ') # just zeros
        else:
            if ( len( Pos.r.offsets) > 2 and equally_spaced( Pos.r.offsets ) ):
                f.write( '\r\n    '+'{:6f}'.format(Pos.r.offsets[1])+'  ')
                f.write( '\r\n    '+'{:6f}'.format(Pos.r.offsets[-1])+'  ')
            elif (len(Pos.r.offsets) == 1):
                f.write( '\r\n    '+'{:6f}'.format(Pos.r.offsets[0])+'  ')
            else:
                for i in range(len(Pos.r.offsets)):
                    f.write( '\r\n    '+'{:6f}'.format(Pos.r.offsets[i]))
        f.write( '/ \t \t ! RR(1)  ... (m) \r\n' )
    return

def fileparts(fname):
    fpath = os.path.dirname(os.path.abspath(fname))
    if '.' not in fname:
        ext = None
    else:
        ext_ind = [i for i in range(len(fname)) if fname[i]=='.']
        if len(ext_ind) > 1:
            raise ValueError("something fishy about a filename with two periods")
        else:
            ext_ind = ext_ind[0]
        ext = fname[ext_ind:]
        fname = fname[:ext_ind]
    return fpath, fname, ext

class Empty():
    def __init__(self):
        return

def write_env( envfil, model, TitleEnv, freq, ssp, bdry, pos, beam, cint, RMax, *varargin ):
    if (envfil[-4:] != '.env' ):
        envfil = envfil + '.env' # append extension

    if ( len( varargin ) == 0 ):
        f = open( envfil, 'w' ) # create new env file
    else:
        f = open( envfil, 'a' )   # append to existing envfil


    model =  model.upper()   # convert to uppercase


    f.write('\'' + TitleEnv + '\' ! Title \r\n')
    f.write('{:8.2f}'.format(freq) +' \t \t \t ! Frequency (Hz) \r\n')
    f.write('{:5d}'.format(ssp.NMedia)+ ' \t \t \t ! NMedia \r\n')
    f.write('\'' + bdry.Top.Opt + '\''+ ' \t \t \t ! Top Option \r\n')

    if ( bdry.Top.Opt[0] == 'A' ): # analytic boundary
        f.write('     {:6.2f}'.format(ssp.depth[0]) + \
                    ' {:6.2f}'.format(bdry.Top.hs.alphaR) + \
                    ' {:6.2f}'.format(bdry.Top.hs.betaR) + \
                    ' {:6.2g}'.format(bdry.Top.hs.rho) + \
                    ' {:6.2f}'.format(bdry.Top.hs.alphaI) + \
                    ' {:6.2f}'.format(bdry.Top.hs.betaI) + \
                    '  \t ! upper halfspace \r\n')


    # SSP
    for medium in range(ssp.NMedia):
        f.write('{:5d}'.format(ssp.N[ medium ]) + \
                ' {:4.2f}'.format(ssp.sigma[ medium ]) + \
                ' {:6.2f}'.format(ssp.depth[ medium+1 ]) + ' \t ! N sigma depth \r\n') 
        for ii in range(len( ssp.raw[ medium].z )):
            f.write('\t {:6.2f} '.format(ssp.raw[ medium ].z[ ii ]) + \
                '{:6.2f} '.format(ssp.raw[ medium ].alphaR[ ii ]) + \
                '{:6.2f} '.format(ssp.raw[ medium ].betaR[ ii ]) + \
                '{:6.2g}'.format(ssp.raw[ medium ].rho[ ii ] ) +  \
                ' {:10.6f} '.format(ssp.raw[ medium ].alphaI[ ii ]) + \
                '{:6.2f} '.format(ssp.raw[ medium ].betaI[ ii ]) + \
                '/ \t ! z c cs rho \r\n')

    # lower halfspace
    f.write('\''+bdry.Bot.Opt + '\'' + ' {:6.2f}'.format(ssp.sigma[1]) + '  \t \t ! Bottom Option, sigma\r\n') # ssp.sigma( 2 ) )

    fmtr = {'all': lambda x: ' {:6.2f}'.format(float(x))}#, 'void': lambda x: ''}
    a = np.array([])
    d = np.array([float(ssp.depth[ssp.NMedia])])
    vls = [d, bdry.Bot.hs.alphaR, bdry.Bot.hs.betaR, bdry.Bot.hs.rho, \
                    bdry.Bot.hs.alphaI, bdry.Bot.hs.betaI]
    strings = [np.array2string(x, formatter = fmtr) for x in vls]
    s = '  {:6.2f} '.format(ssp.depth[ssp.NMedia])
    for string in strings[1:]:
        if string == '':
            string = ' 0.00'
        s += string
    if ( bdry.Bot.Opt[0] == 'A' ):
#        print('{:6.2f}'.format(bdry.Bot.hs.betaI[0]))
        f.write('   ' + s + '  \t  / \t ! lower halfspace \r\n')

    if  model in  [ 'SCOOTER', 'KRAKEN', 'KRAKENC', 'SPARC' ]:
        f.write('{:6.0f} '.format(cint.Low) + '{:6.0f} \t \t ! cLow cHigh (m/s) \r\n'.format(cint.High))   # phase speed limits
        f.write('{:8.2f} \t \t \t ! RMax (km) \r\n'.format(RMax ))    # maximum range

    # source depths
    f.write('{:5d} \t \t \t \t ! NSD'.format(len( pos.s.depth ) ))

    if ( len( pos.s.depth ) >= 2) and equally_spaced( pos.s.depth ):
        f.write('\r\n    {:6f} '.format(pos.s.depth[0]) +' {:6f} '.format(pos.s.depth[-1]))
    else:
        f.write('\r\n    {:6f}  '.format( pos.s.depth[0] ))

    f.write('/ \t ! SD(1)  ... (m) \r\n' )

    # receiver depths

    f.write('{:5d} \t \t \t \t ! NRD'.format(len( pos.r.depth ) ))

    if ( len( pos.r.depth ) >= 2) and equally_spaced( pos.r.depth ) :
        f.write('\r\n    {:6f} '.format(pos.r.depth[0]) + '{:6f} '.format( pos.r.depth[-1]))
    else:
        f.write('\r\n    ')
        for tmp_depth in pos.r.depth:
            f.write('{:6f} '.format(tmp_depth) )

    f.write('/ \t ! RD(1)  ... (m) \r\n' )

    # receiver ranges
    if (model ==  'BELLHOP' ) or  (model == 'FirePE' ) :
        # receiver ranges
        f.write('{:5d} \t \t \t \t ! NRR'.format(len( pos.r.range ) ))
        
        if ( len( pos.r.range ) >= 2) and equally_spaced( pos.r.range ) :
            f.write('\r\n    {:6f} '.format(pos.r.range[0]) + '{:6f} '.format(pos.r.range[-1]))
        else:
            f.write('\r\n    {:6f}  '.format(pos.r.range[0] ))
        f.write('/ \t ! RR(1)  ... (km) \r\n' )
        write_bell(f, beam )

    f.close()

def write_ssp(sspfile, cw, r_arr):
    """
    Input 
    sspfile - string    
        path of file to write
    cw - numpy 2d array
        ssp profiles, columns are the depth profiles for given range
    r_arr - numpy 1d array
        ranges (in km) for updating ssp
    """
    if sspfile[-3:] != 'ssp':
        sspfile += '.ssp'
    with open(sspfile, 'w') as f:
        f.write(str(r_arr.size)+'\r\n')
        for val in r_arr:
            f.write(str(val) + '\t')
        f.write('\r\n')
        for row in range(cw.shape[0]):
            for val in cw[row,:]:
                f.write(str(val) + '\t')
            f.write('\r\n')
        return
    

def write_bathy( btyfile, range_depth_array):
    """
    Input 
    btyfile - string    
        path of file to write
    range_depth_array - numpy 2d array
        first column is ranges (in km), second column is depths
    """
    if btyfile[-3:] != 'bty':
        btyfile += '.bty'
    with open(btyfile, 'w') as f:
        f.write('L\r\n') 
        f.write(str(range_depth_array.shape[0])+'\r\n')
        for row in range(range_depth_array.shape[0]):
            vals = range_depth_array[row,:]
            f.write('  ' + str(vals[0]) + ' ' + str(vals[1]))
    return

def write_bell(f, beam):
    f.write('\''+beam.RunType+'\''+ ' \t \t \t \t ! Run Type \r\n');

    if (beam.Ibeam != None):
        f.write('{:<i}'.format(beam.Nbeams)+'{:<i}'.format(beam.Ibeam) +' \t \t \t \t ! Nbeams Ibeam \r\n');
    else:
        # if this is a ray trace run and the field beam.Nrays exists to use
        # fewer rays in the trace, then use that
        if ( (beam).RunType[0 ] == 'R') and (beam.Nrays != None):
            f.write('{:d}'.format(beam.Nrays) +' \t \t \t \t \t ! Nbeams \r\n');
        else:
            f.write('{:d}'.format(beam.Nbeams )+'\t \t \t \t \t ! Nbeams \r\n')


    f.write('{:f}'.format( beam.alpha[ 0 ])+' {:f} '.format(beam.alpha[ -1 ])+'/ \t \t ! angles (degrees) \r\n');
    f.write('{:f}'.format(beam.deltas)+' {:f}'.format(beam.box.z)+' {:f}'.format(beam.box.r)+'\t ! deltas (m) box.z (m) box.r (km) \r\n');

# Cerveny-style Gaussian beams
    if ( len( beam.RunType ) > 1) and ( 'GBS' not in Beam.RunType[1:1] ) :
        f.write('\''+beam.Type[0:1]+'\''+' {:f}'.format(Beam.epmult)+' {:f}'.format(Beam.rLoop )+' \t \t ! ''Min/Fill/Cer, Sin/Doub/Zero'' Epsmult RLoop (km) \r\n')
        f.write('{:d}'.format(beam.Nimage)+' {:d}'.format(Beam.Ibwin)+'  \t \t \t \t ! Nimage Ibwin \r\n')

def read_shd_bin(*varargin):
    s = Source(0)
    r = Dom(0,0)
    pos = Pos(s,r)
    '''
    Read TL surfaces from a binary Bellhop/Kraken .SHD file
    without having to convert to ASCII first.
    Useage:
    ... = read_shd_bin( filename, xs, ys )
    where (xs, ys) is the source coordinate in km
    (xs, ys) are optional
    Output is a 4-D pressure field p( Ntheta, Nsd, Nrd, Nrr )
    '''
    if (len(varargin) < 1) or (len(varargin) > 3):
        raise ValueError("Can only pass one to three arguments: filename; (filename, xs, ys); or (filename, freq)")

    filename = varargin[0]

    # optional frequency
    if len(varargin) == 2:
        freq = varargin[1]
    else:
        freq = np.NaN

    # optional source (x,y) coordinate
    if len(varargin) >= 3:
       xs = varargin[1]
       ys = varargin[2]
    else:
       xs = np.NaN
       ys = np.NaN

    ##
    f = open( filename, 'rb' )
    if ( f.fileno == -1 ):
        error( 'read_shd_bin.m: No shade file with that name exists' )

    recl     = unpack('<I', f.read(4))[0];     #record length in bytes will be 4*recl
    title    = unpack('80s', f.read(80))

    f.seek(4 * recl); #reposition to end of first record
    PlotType = unpack('10s', f.read(10))
    PlotType = PlotType

    f.seek(2 * 4 * recl); #reposition to end of second record
    Nfreq  = unpack('<I', f.read(4))[0]
    Ntheta = unpack('<I', f.read(4))[0]
    Nsx    = unpack('<I', f.read(4))[0]
    Nsy    = unpack('<I', f.read(4))[0]
    Nsd    = unpack('<I', f.read(4))[0]
    Nrd    = unpack('<I', f.read(4))[0]
    Nrr    = unpack('<I', f.read(4))[0]
    atten  = unpack('<I', f.read(4))[0]
    f.seek(3 * 4 * recl); #reposition to end of record 3
    freqVec = unpack(str(Nfreq) +'d', f.read(Nfreq*8))

    f.seek(4 * 4 * recl) ; #reposition to end of record 4
    pos.theta   = unpack(str(Ntheta) +'f', f.read(4*Ntheta))[0]

    if ( PlotType[ 1 : 2 ] != 'TL' ):
        f.seek(5 * 4 * recl); #reposition to end of record 5
        pos.s.x     = unpack(str(Nsx)+'f',  f.read(Nsx*4))
        f.seek( 6 * 4 * recl); #reposition to end of record 6
        pos.s.y     = unpack(str(Nsy) + 'f', f.read(Nsy*4))
    else:   # compressed format for TL from FIELD3D
        f.seek(5 * 4 * recl, -1 ); #reposition to end of record 5
        pos.s.x     = f.read(2,    'float32' )
        pos.s.x     = linspace( pos.s.x( 1 ), pos.s.x( end ), Nsx )
        
        f.seek(6 * 4 * recl, -1 ); #reposition to end of record 6
        pos.s.y     = f.read(2,    'float32' )
        pos.s.y     = linspace( pos.s.y( 1 ), pos.s.y( end ), Nsy )

    f.seek(7 * 4 * recl); #reposition to end of record 7
    pos.s.depth = unpack(str(Nsd)+'f', f.read(Nsd*4))
    pos.s.depth = np.array(pos.s.depth)

    f.seek(8 * 4 * recl); #reposition to end of record 8
    pos.r.depth = unpack(str(Nrd) + 'f', f.read(Nrd*4))
    pos.r.depth = np.array(pos.r.depth)

    f.seek(9 * 4 * recl); #reposition to end of record 9
    pos.r.range = unpack(str(Nrr) + 'f',f.read(Nrr*4))
    # pos.r.range = pos.r.range';   # make it a row vector
    pos.r.range = np.array(pos.r.range)

    ##
    # Each record holds data from one source depth/receiver depth pair

    if PlotType == 'rectilin  ':
        pressure = np.zeros(( Ntheta, Nsd, Nrd, Nrr ), dtype=np.complex128)
        Nrcvrs_per_range = Nrd
    if PlotType == 'irregular ':
        pressure = np.zeros(( Ntheta, Nsd,   1, Nrr ), dtype=np.complex128)
        Nrcvrs_per_range = 1
    else:
        pressure = np.zeros(( Ntheta, Nsd, Nrd, Nrr ), dtype=np.complex128)
        Nrcvrs_per_range = Nrd

    ##
    if np.isnan( xs ):    # Just read the first xs, ys, but all theta, sd, and rd
        # get the index of the frequency if one was selected
        ifreq = 0
        if not np.isnan(freq):
           freqdiff = [abs( x - freq ) for x in freqVec]
           ifreq = min( freqdiff )

        for itheta in range (Ntheta):
            for isd in range(Nsd):
                # disp( [ 'Reading data for source at depth ' num2str( isd ) ' of ' num2str( Nsd ) ] )
                for ird in range( Nrcvrs_per_range):
                    recnum = 10 + ( ifreq   ) * Ntheta * Nsd * Nrcvrs_per_range + \
                                  ( itheta  )          * Nsd * Nrcvrs_per_range + \
                                  ( isd     )                * Nrcvrs_per_range + \
                                    ird    
                    status = f.seek(int(recnum) * 4 * recl); #Move to end of previous record
                    if ( status == -1 ):
                        raise ValueError( 'Seek to specified record failed in read_shd_bin' )
                    temp = unpack(str(2*Nrr)+'f', f.read(2 * Nrr*4));    #Read complex data
                    pressure[ itheta, isd, ird, : ] = temp[ 0 : 2 * Nrr -1 : 2 ] + complex(0,1) *np.array((temp[ 1 : 2 * Nrr :2]))
                    # Transmission loss matrix indexed by  theta x sd x rd x rr
                    
    else:              # read for a source at the desired x, y, z.
        
        xdiff = abs( pos.s.x - xs * 1000. )
        [ holder, idxX ] = min( xdiff )
        ydiff = abs( pos.s.y - ys * 1000. )
        [ holder, idxY ] = min( ydiff )
        
        # show the source x, y that was found to be closest
        # [ pos.s.x( idxX ) pos.s.y( idxY ) ]
        for itheta in range(Ntheta):
            for isd in range(Nsd):
                # disp( [ 'Reading data for source at depth ' num2str( isd ) ' of ' num2str( Nsd ) ] )
                for ird in range(Nrcvrs_per_range):
                    recnum = 10 + ( idxX   - 1 ) * Nsy * Ntheta * Nsd * Nrcvrs_per_range +   \
                                  ( idxY   - 1 )       * Ntheta * Nsd * Nrcvrs_per_range +  \
                                  ( itheta - 1 )                * Nsd * Nrcvrs_per_range +  \
                                  ( isd    - 1 )                      * Nrcvrs_per_range + ird - 1
                    status = f.seek(recnum * 4 * recl); # Move to end of previous record
                    if ( status == -1 ):
                        raise ValueError( 'Seek to specified record failed in read_shd_bin' )
                    
                    temp = f.read(2 * Nrr, 'float32' );    #Read complex data
                    pressure[ itheta, isd, ird, : ] = temp[ 1 : 2 : 2 * Nrr ] + complex(0,1) * np.array(temp[ 2 : 2 : 2 * Nrr ])
                    # Transmission loss matrix indexed by  theta x sd x rd x rr
                    
    f.close()
    return [ title, PlotType, freqVec, atten, pos, pressure ] 

def read_shd (*varargin ):
    '''
     Read the shade file
    [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] return vals
     calls the appropriate routine (binary, ascii, or mat file) to read in the pressure field
    
     usage: [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] = read_shd( filename )
        Reads first source.
            [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] = read_shd( filename, xs, ys )
        Reads source at the specified xs, ys coordinate.
            [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] = read_shd( filename, freq )
        Reads source at the specified frequency.
    
     Recommended to include a file extension, if it exists.
     Otherwise it may find a different file than you intended.
    
     Output is a 5-D pressure field p( Nfreq, Ntheta, Nsd, Nrd, Nrr )
    
     If omitted, take a guess at the extension
     Matlab 'exist' command is simpler; however, it searches the whole Matlab search path.
    '''


    # Determine type of file:

    #error( nargchk( 1, 3, len(varargin), 'struct' ) )
    if (len(varargin) < 1) or (len(varargin) > 3):
        raise ValueError("Can only pass one to three arguments: filename; (filename, xs, ys); or (filename, freq)")

    filename = varargin[0]

    # optional frequency
    if len(varargin) == 2:
        freq = varargin[1]
    else:
        freq = np.NaN

    # optional source (x,y) coordinate
    if len(varargin) >= 3:
       xs = varargin[1]
       ys = varargin[2]
    else:
       xs = np.NaN
       ys = np.NaN

    PlotType = [];  # in case this was not set

    [holder , name, ext ] = fileparts( filename )
    if (  ext == '.mat' ) :
        [ holder, holder1, ext2 ] = fileparts( name )

        if ext2 == '.shd':
            FileType = 'shdmat'
        elif ext2 == '.grn':
            FileType = 'grnmat'
        else:
            pass
    else:
        if filename == 'ASCFIL':
             FileType = 'asc'
        elif filename == 'SHDFIL':
            FileType = 'shd'
        elif filename == 'tl.grid':
            FileType = 'RAM'
        else:
            endchar = len( filename )
            if ( endchar >= 4 ):
                FileType = filename[ endchar - 3 : endchar ].lower()

    ##
    if FileType in ['shd', 'grn' ]:   # binary format
        if len(varargin) ==  1:
            [ PlotTitle, PlotType, freqVec, atten, pos, pressure ] = read_shd_bin( filename )
        if len(varargin) ==  2:
            [ PlotTitle, PlotType, freqVec, atten, pos, pressure ] = read_shd_bin( filename, freq )
        if len(varargin) ==  3:
            [ PlotTitle, PlotType, freqVec, atten, pos, pressure ] = read_shd_bin( filename, xs, ys )
    elif FileType ==  'shdmat':  # Shade function mat file
        loadmat( filename )

        # has a specific source xs, ys been given?
        if not np.isnan( xs ):
            xdiff = abs( pos.s.x - xs * 1000. )
            [ holder, idxX ] = min( xdiff )
            ydiff = abs( pos.s.y - ys * 1000. )
            [ holder, idxY ] = min( ydiff )

            # extract the appropriate source index
            pressureT = pressure[ idxX, idxY, :, :, :, : ]
            pressure = np.reshape( pressureT, [ pos.Nsd, pos.Nrd, pos.Ntheta, pos.Nrr ] )

        # has a specific frequency been given?
        if not np.isnan( freq ):
            freqdiff = abs( freqVec - freq )
            [ holder, ifreq ] = min( freqdiff )
            # extract the appropriate source index
            pressure = pressure[ ifreq, 1, :, : ]

    elif FileType ==  'asc': # ascii format
        [ PlotTitle, PlotType, freqVec, atten, pos, pressure ] = read_shd_asc( filename )

    elif FileType ==  'grnmat':   # Green's function mat file
        loadmat( filename )
        pos.r.range = np.array(pos.r.range.T);   # make it a column vector to match read_shd_bin

    elif FileType ==  'RAM':
        [ PlotTitle, PlotType, freqVec, atten, pos, pressure ] = read_ram_tlgrid
    else:
        raise ValueError( 'Unrecognized file extension' )

    # clean up PlotTitle by taking only the part up inside the quotes
    # nchars = strfind( PlotTitle, '''' );   # find quotes
    # PlotTitle = [ PlotTitle( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) ]
    return [ PlotTitle, PlotType, freqVec, atten, pos, pressure ] 

def crci( c, alpha, freq, AttenUnit ):
    '''
     Convert real wave speed and attenuation to a single complex wave speed
     5 cases:
     N for Nepers/meter
     M for dB/meter      (M for Meters)
     F for dB/m-kHZ      (F for frequency dependent)
     W for dB/wavelength (W for Wavelength)
     Q for Q
     T for Thorpe
    '''
    omega = 2.0 * np.pi * freq

    # *** Convert to Nepers/m ***

    alphaT = 0.0

    if AttenUnit[0] == 'N' :   # Nepers/m
       alphaT = alpha
    elif AttenUnit[0] == 'M' :   # dB/meter
       alphaT = alpha / 8.6858896
    elif AttenUnit[0] == 'F' :   # dB/m-kHZ
       alphaT = alpha * freq / 8685.8896
    elif AttenUnit[0] == 'W' :   # dB/wavelength
       if ( c != 0.0 ):
           alphaT = alpha * freq / ( 8.6858896 * c )
    elif AttenUnit[0] == 'Q' :
       if( c * alpha != 0.0 ):
           alphaT = omega / ( 2.0 * c * alpha )
    else:
        print('Unknown attenuation unit')

    # added volume attenuation
    if AttenUnit[1] ==  ( 'T' ) :  # Thorp
        f2     = np.square(( freq / 1000.0 ))
        # Original Thorp (1967) formula
        #    alphaT = 40.0 * f2 / ( 4100.0 + f2 ) + 0.1 * f2 / ( 1.0 + f2 );
        #    alphaT = alphaT / 914.4;     # dB / m
        #    alphaT = alphaT / 8.6858896; # Nepers / m
        
        # Updated formula from JKPS Eq. 1.34
        Thorp = 3.3e-3 + 0.11 * f2 / ( 1.0 + f2 ) + 44.0 * f2 / ( 4100.0 + f2 ) + 3e-4 * f2  # dB/km
        Thorp = Thorp / 8685.8896 # Nepers / m
        alphaT = alphaT + Thorp
    #     *** Convert Nepers/m to equivalent imaginary sound speed ***
    alphaT = alphaT * c*c / omega
    crci   = complex(c, alphaT)
    return crci

def topbot( lines, line_ind, freq, BCType, AttenUnit ):
    ''' Handles top and bottom boundary conditions
    Input:
    ENVFIL: Environmental file
    freq:   frequency    #   BCType: Boundary condition type

    Output:
    Bdry.cp:    P-wave speed in halfspace
    Bdry.cs:    S-wave speed in halfspace
    Bdry.rho:   density in halfspace

    BumDen:  Bump density
    eta:     Principal radius 1
    xi:      Principal radius 2
    '''
    global alphaR, betaR, rhoR, alphaI, betaI

    # *** Echo to PRTFIL user's choice of boundary condition ***
    bctypes = ['S', 'H', 'T', 'I', 'V', 'R', 'A', 'F', 'W', 'P']
    if BCType not in bctypes:
       raise ValueError( 'Fatal error: Unknown boundary condition type' )

    # ****** Read in BC parameters depending on particular choice ******

    cp  = 0.0
    cs  = 0.0
    rho = 0.0
    ice = Ice(0.0, 0.0, 0.0)
    hs = None

    # *** Twersky ice model parameters ***

    if ( BCType == 'S' or BCType == 'H' or BCType == 'T' or BCType == 'I' ):
        tmp = lines[line_ind]
        ice.BumDen, ice.eta, ice.xi = float(tmp.split()[0]), float(tmp.split()[1]), float(tmp.split()[2])
        line_ind += 1

    # *** Half-space properties ***

    if ( BCType == 'A' ):
        tmp = lines[line_ind]
        line_ind += 1
        tmp = tmp.split()
        tmp = [x for x in tmp if x[0].isnumeric()] # filter out non numbers
        num_vals = len(tmp)
        if num_vals == 6:             
            ztmp, alphaR, betaR, rhoR, alphaI, betaI = [float(x) for x in tmp[0:6]]
        elif num_vals == 5:
            ztmp, alphaR, betaR, rhoR, alphaI = [float(x) for x in tmp[0:5]]
        elif num_vals == 4:
            ztmp, alphaR, betaR, rhoR = [float(x) for x in tmp[0:4]]
        elif num_vals == 3:
            ztmp, alphaR, betaR = [float(x) for x in tmp[0:3]]
        elif num_vals == 2:
            ztmp, alphaR = [float(x) for x in tmp[0:2]]
        elif num_vals == 1:
            ztmp = [float(x) for x in tmp][0]
        else: # there were no vals to read in so defaults will be used
            pass

        cp  = crci( alphaR, alphaI, freq, AttenUnit )
        cs  = crci( betaR,  betaI,  freq, AttenUnit )
        rho = rhoR
       
        # raw values in user units
        hs = HS(alphaR, betaR, rho, alphaI, betaI)
    return [ cp, cs, rho, ice, hs, line_ind]

def read_env_core( envfil ):
    global HV, NFirstAcoustic, NLastAcoustic
    global alphaR, betaR, rhoR, alphaI, betaI

    alphaR = 1500   # defaults
    betaR  = 0
    rhoR   = 1
    alphaI = 0
    betaI  = 0

    bdry = Bndry(None,None)
    bdry.Top = TopBndry('dummy_opt')
    bdry.Bot = BotBndry('dummy_opt', 10)

    NFirstAcoustic = 0
    bdry.Top.cp    = 0.0
    bdry.Top.cs    = 0.0
    bdry.Top.rho   = 0.0
    bdry.Bot.cp    = 2000.0
    bdry.Bot.cs    = 0.0
    bdry.Bot.rho   = 2.0

    ssp = None

    with open(envfil, 'r') as f:
        lines = f.readlines()
        # first line is title
        line_ind = 0
        TitleEnv = re.findall('\'(.*)\'',lines[line_ind] )[0]
        # second line is freq
        line_ind += 1
        freq     = float((lines[line_ind].split())[0])
        line_ind += 1

        # third line is number of media
        tmp = lines[line_ind]
        line_ind += 1
        NMedia   = int(tmp.split()[0])
        # initialize depth array
        ssp_depth = [0]*(NMedia+1)
        
        # fourth line is Top option 
        tmp = lines[line_ind]    
        line_ind += 1
        TopOpt   = re.findall('\'(.*)\'', tmp)[0]
        if len(TopOpt) < 5:
            TopOpt = TopOpt + ' '*(4-len(TopOpt)) # whitespace pad to 4 characters
        # convert the deprecated '*' option to '~'
        bdry.Top.Opt = TopOpt.replace('*', '~')

        SSPType = bdry.Top.Opt[0]
        bdry.Top.BC     = bdry.Top.Opt[1]
        # this should throw error, I may need to pad the end with spaces
        AttenUnit = bdry.Top.Opt[2:4]

        #     *** SSP approximation options ***
        ssp_types = ['N', 'C', 'S', 'Q', 'H', 'A']
        if SSPType not in ssp_types:
            raise ValueError('SSP must be in ' + str(ssp_types))


        #     *** Attenuation options ***
        atten_opts = ['N', 'F', 'M', 'W', 'Q']
        if AttenUnit[0] not in atten_opts:
            raise ValueError( 'Unknown attenuation units, must be ' + str(atten_opts) )
        
        [ bdry.Top.cp, bdry.Top.cs, bdry.Top.rho, holder, bdry.Top.hs, line_ind ] = topbot(lines, line_ind, freq, bdry.Top.BC, AttenUnit )

        ssp_z   = []
        ssp_c   = []
        ssp_cs  = []
        ssp_rho = []
        ssp_raw_list = []
        Loc = [0] * NMedia
        ssp_Npts = [0]*NMedia
        ssp_N = [0]*NMedia
        ssp_sigma = [0]*(NMedia+1)
        for medium in range(NMedia):
            if ( medium == 1 ):
                Loc[medium] = 0
            else:
                Loc[ medium ] = Loc[ medium - 1 ] + ssp_Npts[ medium - 1 ]
               
            # next line is the num points in ssp, sigma, depth
            tmp = lines[line_ind]
            line_ind += 1
            ssp_N[medium], ssp_sigma[medium], ssp_depth[medium+1] = int(float(tmp.split()[0])), float(tmp.split()[1]), float(tmp.split()[2])
           
            # read in the SSP
            ssp_raw_list.append(SSPraw([], [], [], [], [], []))
           
            for ii in range(1000):
                tmp = lines[line_ind]
                line_ind += 1
                tmp = tmp.split()
                tmp = [x for x in tmp if x[0].isnumeric()] # filter out non numbers
                num_vals = len(tmp)
                if num_vals == 6:             
                    ztmp, alphaR, betaR, rhoR, alphaI, betaI = [float(x) for x in tmp[0:6]]
                elif num_vals == 5:
                    ztmp, alphaR, betaR, rhoR, alphaI = [float(x) for x in tmp[0:5]]
                elif num_vals == 4:
                    ztmp, alphaR, betaR, rhoR = [float(x) for x in tmp[0:4]]
                elif num_vals == 3:
                    ztmp, alphaR, betaR = [float(x) for x in tmp[0:3]]
                elif num_vals == 2:
                    ztmp, alphaR = [float(x) for x in tmp[0:2]]
                elif num_vals == 1:
                    ztmp = [float(x) for x in tmp[0]][0]
                else: # there were no vals to read in so defaults will be used
                    pass

             
                ## TRANSLATE CRCI, IN BELLHOP DIR 
                cp = crci( alphaR, alphaI, freq, AttenUnit )
                cs = crci( betaR,  betaI,  freq, AttenUnit )
                ssp_z.append(ztmp)   # add in to existing vector
                ssp_c.append(cp)
                ssp_cs.append(cs)
                ssp_rho.append(rhoR)

                ssp_raw_list[ medium].z.append(ztmp)   # add in to existing vector
                ssp_raw_list[ medium].alphaR.append(alphaR)
                ssp_raw_list[ medium].alphaI.append(alphaI)
                ssp_raw_list[ medium].betaR.append(betaR)
                ssp_raw_list[ medium].betaI.append(betaI)
                ssp_raw_list[ medium].rho.append(rhoR)

                # check for end of this layer
                if ( ztmp == ssp_depth[ medium+1 ] ):
                    ssp_depth[0] = ssp_z[0]
                    # break
                    break
         
            # end of loop for that medium 
            if ( ssp_N[ medium ] == 0 ): # calculate mesh automatically
                # choose a reference sound speed
                C = alphaR
                if ( betaR > 0.0 ):
                    C = betaR  # shear?
                deltaz = 0.05 * C / freq     # default sampling: 20 points per wavelength
                ssp_N[ medium ] = round( ( ssp_depth[ medium + 1 ] - ssp_depth[ medium ] ) / deltaz )
                ssp_N[ medium ] = np.max( [ssp_N[ medium ], 10] )     # require a minimum of 10 points
           
            # keep track of first and last acoustic medium
            if ( not np.any( cs ) ):   # shear anywhere?
                if ( NFirstAcoustic == 0 ):
                    NFirstAcoustic = medium
                NLastAcoustic  = medium
       
            # stuff for Bellhop
            if ( medium == 0 ):
                HV       = np.diff( np.array(ssp_z) ) # layer thicknesses
                ssp_cz   = np.diff( np.array(ssp_c )) / HV # gradient of ssp (centered-difference approximation)
            ssp_Npts[ medium ] = ii
            if ( medium == 0 ):
                ssp_depth[0] = ssp_z[0]

        ## lower halfspace

        # read in next line
        tmp = lines[line_ind]
        line_ind += 1
        BotOpt   = re.findall('\'(.*)\'', tmp)[0]

        if len(BotOpt) < 5:
            BotOpt = BotOpt + ' '*(4-len(BotOpt)) # whitespace pad to 4 characters
        # convert the deprecated '*' option to '~'
        bdry.Bot.Opt = BotOpt.replace('*', '~')
        bdry.Bot.BC    = bdry.Bot.Opt[0]

        ## TRANSLATE TOPBOT IN BELLHOP DIR
        [ bdry.Bot.cp, bdry.Bot.cs, bdry.Bot.rho, holder, bdry.Bot.hs, line_ind ] = topbot( lines, line_ind, freq, bdry.Bot.BC, AttenUnit )
            

        bdry.Top.depth = ssp_depth[0]
        bdry.Bot.depth = ssp_depth[ NMedia ]

        # Get rho, c just INSide the boundary (used later for reflection
        # coefficients)
        I = NFirstAcoustic
        bdry.Top.rhoIns = ssp_rho[ I ]
        bdry.Top.cIns   = ssp_c[   I ]

        I = Loc[ NLastAcoustic ] + ssp_Npts[ NLastAcoustic ]
        bdry.Bot.rhoIns = ssp_rho[ I ]
        bdry.Bot.cIns   = ssp_c[   I ]
        ssp = SSP(ssp_raw_list, ssp_depth, NMedia, N=ssp_N, sigma=ssp_sigma)
    return [TitleEnv, freq, ssp, bdry, lines, line_ind]

def readvector(lines, line_ind):
    '''
    Reads a row-vector

    This routine emulates a Fortran language capability that allows a '/'
    to be used to terminate input

    A user can type:
      5     ! Number of values
      0 1000 /
    as a shorthand for
      5
      0 250 500 750 1000
    Here the '/' terminates the reading and the code then creates 4 equally
    spaced values between 0 and 1000.
    However, one also has the option of using specific (non-equispaced
    values)
    '''
    line_ind = int(line_ind)
    tmp = lines[line_ind]
    line_ind += 1
    Nx = float(tmp.split()[0])
    tmp = lines[line_ind]
    line_ind += 1
    
    if '/' in tmp:
        tmp = tmp.split()
        tmp = [x for x in tmp if x[0].isnumeric()]
        if len(tmp) == 1:
            x = float(tmp[0])
        else:
            x = [float(tmp[0]), float(tmp[1])]
    if Nx > 2:
        x = np.linspace( x[ 0 ], x[ 1 ], Nx ) # generate the vector
    else:
        x = [float(x) for x in tmp if x[0].isnumeric()]
    return x, Nx, line_ind

def readsdrd(lines, line_ind):
    '''
     Read source depths and receiver depths

     Variable 'Pos' is a structure:
     Pos.r.depth = vector of receiver depths
     Pos.Nrd     = number of receiver depths
     Pos.s.depth = vector of source depths
     Pos.Nsd     = number of source depths
    '''
    [ depth, Nx, line_ind] = readvector(lines, line_ind);
    source = Source(depth)

    [ depth, Nrd, line_ind ] = readvector(lines, line_ind);
    dom = Dom(None, depth)
    pos = Pos(source, dom) 
    return pos, line_ind

def  read_env( envfil, model ):
    ''' Read an environmental file

     should check here that model is one of the known ones
    '''
    if (  envfil != 'ENVFIL' ) and  (envfil[-4:] != '.env' ):
        envfil =  envfil + '.env' # append extension

    model = model.upper()    # convert to uppercase

    [ TitleEnv, freq, ssp, bdry, lines, line_ind] = read_env_core( envfil )     # read in the environmental file
       
    model_cat =  ['SCOOTER', 'KRAKEN', 'KRAKENC', 'KRAKEL', 'SPARC' ]

    if model in model_cat:
        tmp = lines[line_ind]
        line_ind += 1
        tmp = [x for x in tmp.split() if x[0].isnumeric()]
        low, high = [float(x) for x in tmp]
        cint = cInt(low, high)
        tmp = lines[line_ind]
        line_ind += 1 
        RMax = float(tmp.split()[0])

        
    else:   # dummy values for BELLHOP
        cint = cInt(1500, 1e9)
        RMax      = 100 

    # BELLHOP3D has x-y coordinates of sources as well

   
    if ( model == 'BELLHOP3D' ) :
        raise ValueError("Bellhop 3d not supported atm")

    # !!! check: does this delete Pos.s.x, etc. when executed
    [pos, line_ind] = readsdrd(lines, line_ind)                            # read in the source and receiver depths

    if (  model ==  'BELLHOP' ) :
        raise ValueError("Bellhop not supported at the moment")
#        Pos.r.range = readr( fid )      # read in receiver ranges
#        Pos.Nrr     = length( Pos.r.range ) 
#        Beam        = read_bell( fid, Bdry, freq, Bdry.Bot.depth, Bdry.Top.depth, Pos.r.range( end ) ) 
#    elseif ( strcmp( model, 'BELLHOP3D' ) )
#        Pos.r.range = readr( fid )      # read in receiver ranges
#        Pos.Nrr     = length( Pos.r.range ) 
#        Pos.theta   = readtheta( fid ) 
#        Pos.Ntheta  = length( Pos.theta ) 
#        Beam        = read_bell( fid, Bdry, freq, Bdry.Bot.depth, Bdry.Top.depth, Pos.r.range( end ) ) 
    else:   # dummy value for models that don't use Beam parameters
        beam = Beam('CG', Nbeams=0, alpha= [-15,15], box=Box(1.05*RMax, 1.05*np.max(pos.r.depth)), deltas = 0) 
        pos.r.range     = np.linspace( 0, RMax, 501 )    # set up receiver range vector
    return  [ TitleEnv, freq, ssp, bdry, pos, beam, cint, RMax]

def read_modes(**kwargs):
    '''
     Read the modes produced by KRAKEN
     usage:
     keys are 'fname', 'freq', 'modes'
    'fname' and 'freq' are mandatory, 'modes' is if you only want a subset of modes
        [ Modes ] = read_modes_bin( filename, modes )
     filename is without the extension, which is assumed to be '.moA'
     freq is the frequency (involved for broadband runs)
       (you can use freq=0 if there is only one frequency)
     modes is an optional vector of mode indices

     derived from readKRAKEN.m    Feb 12, 1996 Aaron Thode

     Translated to python by Hunter Akins 2019

     Modes.M          number of modes
     Modes.k          wavenumbers
     Modes.z          sample depths for modes
     Modes.phi        modes

     Modes.Top.bc
     Modes.Top.cp
     Modes.Top.cs
     Modes.Top.rho
     Modes.Top.depth

     Modes.Bot.bc
     Modes.Bot.cp
     Modes.Bot.cs
     Modes.Bot.rho
     Modes.Bot.depth

     Modes.N          Number of depth points in each medium
     Modes.Mater      Material type of each medium (acoustic or elastic)
     Modes.Nfreq      Number of frequencies
     Modes.Nmedia     Number of media
     Modes.depth      depths of interfaces
     Modes.rho        densities in each medium
     Modes.freqVec    vector of frequencies for which the modes were calculated
    '''

    
    filename = kwargs['fname']
    if 'freq' in kwargs.keys():
        freq = kwargs['freq']
    if 'modes' in kwargs.keys():
        modes = kwargs['modes']
    with open(filename, 'rb') as f:
        iRecProfile = 1;   # (first time only)
        
        lrecl     = 4*unpack('<I', f.read(4))[0];     #record length in bytes

        rec = iRecProfile - 1;

        f.seek(rec * lrecl + 4) # do I need to do this ?

        title    = unpack('80s', f.read(80))
        Nfreq  = unpack('<I', f.read(4))[0]
        Nmedia = unpack('<I', f.read(4))[0]
        Ntot = unpack('<l', f.read(4))[0]
        Nmat = unpack('<l', f.read(4))[0]
        N = []
        Mater = []


        if Ntot < 0:
            return

        # N and Mater
        rec   = iRecProfile;
        f.seek(rec * lrecl); # reposition to next level
        for Medium in range(Nmedia):
           N.append(unpack('<I', f.read(4))[0])
           Mater.append(unpack('8s', f.read(8))[0])


        # depth and density
        rec = iRecProfile + 1
        f.seek(rec * lrecl)
        bulk        = unpack('f'*2*Nmedia, f.read(4*2*Nmedia))
        depth = [bulk[i] for i in range(0,2*Nmedia,2)]
        rho = [bulk[i+1] for i in range(0,2*Nmedia,2)]

        # frequencies
        rec = iRecProfile + 2;
        f.seek(rec * lrecl);
        freqVec = unpack('d', f.read(8))[0]
        freqVec = np.array(freqVec)

        # z
        rec = iRecProfile + 3
        f.seek(rec * lrecl)
        z = unpack('f'*Ntot, f.read(Ntot*4))

        # read in the modes

        # identify the index of the frequency closest to the user-specified value
        freqdiff = abs(freqVec - freq );
        freq_index = np.argmin( freqdiff );

        # number of modes, m
        iRecProfile = iRecProfile + 4;
        rec = iRecProfile;

        # skip through the mode file to get to the chosen frequency
        for ifreq in range(freq_index+1):
            f.seek(rec * lrecl);
            M = unpack('l', f.read(8))[0]
       
           
           # advance to the next frequency
            if ( ifreq < freq_index ):
                iRecProfile = iRecProfile + 2 + M + 1 + floor( ( 2 * M - 1 ) / lrecl );   # advance to next profile
                rec = iRecProfile;
                f.seek(rec * lrecl)

        if 'modes' not in kwargs.keys():
            modes = np.linspace(0, M-1, M, dtype=int);    # read all modes if the user didn't specify

        # Top and bottom halfspace info

        # Top
        rec = iRecProfile + 1
        f.seek(rec * lrecl)
        top_bc    = unpack('c', f.read(1))[0]
        cp              = unpack('ff', f.read(8))
        top_cp    = complex( cp[ 0 ], cp[ 1 ] )
        cs              = unpack('ff', f.read(8))
        top_cs    = complex( cs[ 1 ], cs[ 1 ] )
        top_rho   = unpack('f', f.read(4))[0]
        top_depth = unpack('f', f.read(4))[0]
    
        top_hs = HS(alphaR=top_cp.real, alphaI=top_cp.imag, betaR=top_cs.real, betaI=top_cs.imag)
        top = TopBndry(top_bc, depth=top_depth)  

        # Bottom
        bot_bc    = unpack('c', f.read(1))[0]
        cp              = unpack('ff', f.read(8))
        bot_cp    = complex( cp[ 0 ], cp[ 1 ] )
        cs              = unpack('ff', f.read(8))
        bot_cs    = complex( cs[ 1 ], cs[ 1 ] )
        bot_rho   = unpack('f', f.read(4))[0]
        bot_depth = unpack('f', f.read(4))[0]
  
        bot_hs = HS(alphaR=bot_cp.real, alphaI=bot_cp.imag, betaR=bot_cs.real, betaI=bot_cs.imag)
        bot = BotBndry(bot_bc, bot_hs, depth=bot_depth)  

        rec = iRecProfile
        f.seek(rec * lrecl)
        # if there are modes, read them
        if ( M == 0 ):
           modes_phi = []
           modes_k   = []
        else:
            modes_phi = np.zeros((Nmat, len( modes )),dtype=np.complex64)   # number of modes
           
            for ii in range(len(modes)):
                rec = iRecProfile + 2 + int(modes[ ii ])
                f.seek(rec * lrecl)
                phi = unpack('f'*2*Nmat, f.read(2*Nmat*4)) # Data is read columwise
                phir = np.array([phi[i] for i in range(0,2*Nmat,2)])
                phii = np.array([phi[i+1] for i in range(0,2*Nmat,2)])
                
                modes_phi[ :, ii ] = phir + complex(0, 1)*phii;
           
            rec = iRecProfile + 2 + M;
            f.seek(rec * lrecl)
            k    = unpack('f'*2*M, f.read(4*2*M))
            kr = np.array([k[i] for i in range(0,2*M,2)])
            ki = np.array([k[i+1] for i in range(0,2*M,2)])
            modes_k = kr+ complex(0,1) * ki
            modes_k = np.array([modes_k[i] for i in modes], dtype=np.complex64)  # take the subset that the user specified

    input_dict = {'M':M, 'modes_k': modes_k, 'z':z, 'modes_phi':modes_phi, 'top':top, 'bot':bot, 'N':N, 'Mater':Mater, 'Nfreq': Nfreq, 'Nmedia': Nmedia, 'depth':depth, 'rho':rho, 'freqVec':freqVec}
    modes = Modes(**input_dict)
    return modes

def my_float(x):
    if x == 'nan' or x == 'NaN':
        return np.nan
    else:
        return float(x)

def read_arrivals_asc(fname, narrmx=200):
#function [ Arr, Pos ] = read_arrivals_asc( ARRFile, Narrmx )

# Read the arrival time/amplitude data computed by BELLHOP
#
# usage:
#[ Arr, Pos ] = read_arrivals_asc( ARRFile, Narrmx );
#
# Arr is a structure containing all the arrivals information
# Pos is a structure containing the positions of source and receivers
#
# ARRFile is the name of the Arrivals File
# Narrmx is the maximum number of arrivals allowed
# mbp 9/96
    if fname[-4:] != '.arr':
        fname = fname + '.arr'
    with open(fname, 'r') as f:
        # read the header info
        lines = f.readlines()
        tmp = lines[0]
        a = tmp.split()
        freq = float(a[0])
        Nsd, Nrd, Nrr = int(a[1]), int(a[2]), int(a[3])
        tmp = lines[1].split()
        sd = [float(x) for x in tmp]
        tmp = lines[2].split()
        rd = [float(x) for x in tmp]
        tmp = lines[3].split()
        rr = [float(x) for x in tmp] # in meters
        dom = Dom(np.array(rr), np.array(rd))
        arrival_list = []
        source = Source(np.array(sd))
        pos = Pos(source, dom)
        line_index= 4
        for i in range(Nsd):
            num_angles = int(lines[line_index].split()[0])
            line_index += 1
            for j in range(Nrd):
                for k in range(Nrr):
                    num_arrivals = int(lines[line_index].split()[0])
                    line_index += 1
                    loc_arrivals = [] #arrivals for this specific sd, rd, and rr
                    for arr in range(num_arrivals):
                        tmp = lines[line_index].split()
                        amp = my_float(tmp[0])*np.exp(complex(0,1)*my_float(tmp[1])*np.pi/180)
                        delay = complex(my_float(tmp[2]), my_float(tmp[3]))
                        src_ang = my_float(tmp[4])
                        rec_ang = my_float(tmp[5])
                        num_top_bnc = int(tmp[6])
                        num_bot_bnc = int(tmp[7])
                        loc_arrivals.append([amp, delay, src_ang, rec_ang, num_top_bnc, num_bot_bnc])
                        line_index += 1
                    arrival_list.append(loc_arrivals)
    return arrival_list, pos
    #% loop to read all the arrival info (delay and amplitude)

    #Arr.A         = zeros( Nrr, Narrmx, Nrd, Nsd );
    #Arr.delay     = zeros( Nrr, Narrmx, Nrd, Nsd );
    #Arr.SrcAngle  = zeros( Nrr, Narrmx, Nrd, Nsd );
    #Arr.RcvrAngle = zeros( Nrr, Narrmx, Nrd, Nsd );
    #Arr.NumTopBnc = zeros( Nrr, Narrmx, Nrd, Nsd );
    #Arr.NumBotBnc = zeros( Nrr, Narrmx, Nrd, Nsd );

    #for isd = 1 : Nsd
    #   Narrmx2 = fscanf( fid, '%i', 1 );  % max. number of arrivals to follow
    #   disp( [ 'Max. number of arrivals for source index ', num2str( isd ), ' is ', num2str( Narrmx2 ) ] );
    #   for ird = 1:Nrd
    #      for ir = 1:Nrr
    #         Narr = fscanf( fid, '%i', 1 );	% number of arrivals
    #         Arr.Narr( ir, ird, isd ) = Narr;
    #         
    #         if Narr > 0   % do we have any arrivals?
    #            da = fscanf( fid, '%f', [ 8, Narr ] );
    #            Narr = min( Narr, Narrmx ); % we'll keep no more than Narrmx values
    #            Arr.Narr( ir, ird, isd ) = Narr;

    #            Arr.A(         ir, 1:Narr, ird, isd ) = da( 1, 1:Narr ) .* exp( 1i * da( 2, 1:Narr ) * pi/180);
    #            Arr.delay(     ir, 1:Narr, ird, isd ) = da( 3, 1:Narr ) + 1i * da( 4, 1:Narr );
    #            Arr.SrcAngle(  ir, 1:Narr, ird, isd ) = da( 5, 1:Narr );
    #            Arr.RcvrAngle( ir, 1:Narr, ird, isd ) = da( 6, 1:Narr );
    #            Arr.NumTopBnc( ir, 1:Narr, ird, isd ) = da( 7, 1:Narr );
    #            Arr.NumBotBnc( ir, 1:Narr, ird, isd ) = da( 8, 1:Narr );
    #         end
    #      end		% next receiver range
    #   end		% next receiver depth
    #end	% next source depth

    #fclose( fid );
