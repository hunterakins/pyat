import numpy as np
from scipy.io import loadmat
import os
from struct import unpack
from pyat import objects



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
def write_fieldflp(flpfil, Option, Pos):
    if flpfil[-4:] != '.flp':
        flpfil = flpfil + '.flp'

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
            f.write('\r\n    '+'{:6f}'.format(Pos.s.depth[ 1 ]+'  '))
            f.write('\r\n    '+'{:6f}'.format(Pos.s.depth[ -1 ]+'  '))
        elif len(Pos.s.depth) == 1:
            f.write( '\r\n    '+ '{:6f}'.format(Pos.s.depth[0])+'  ')
        else:
            for i in range(len(Pos.s.depth)):
                f.write( '\r\n    '+ '{:6f}'.format(Pos.s.depth[i])+'  ')

        f.write( '/ \t ! SD(1)  ... (m) \r\n' )

        # receiver depths

        f.write('{:5}'.format(len(Pos.r.depth)) + ' \t \t \t \t ! NRD')
    
        if ( len( Pos.r.depth ) > 2 and equally_spaced( Pos.r.depth ) ):
            f.write( '\r\n    '+'{:6f}'.format(Pos.r.depth[1])+'  ')
            f.write( '\r\n    '+'{:6f}'.format(Pos.r.depth[-1])+'  ')
        elif (len(Pos.r.depth) == 1):
            f.write( '\r\n    '+'{:6f}'.format(Pos.r.depth[0])+'  ')
        else:
            for i in range(len(Pos.r.depth)):
                f.write( '\r\n    '+'{:6f}'.format(Pos.r.depth[i]))

        f.write( '/ \t ! RD(1)  ... (m) \r\n' )

        # receiver range offsets

        f.write( '{:5}'.format(len( Pos.r.depth ))+' \t \t ! NRR',  )
        f.write( '\r\n    '+'{:6.2f}'.format(0.0) +'  ')
        f.write( '\r\n    '+'{:6.2f}'.format(0.0) +'  ')
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

def read_shd_bin(*varargin):
    Pos = Empty()
    Pos.s = Empty()
    Pos.r = Empty()

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
    Pos.theta   = unpack(str(Ntheta) +'f', f.read(4*Ntheta))[0]

    if ( PlotType[ 1 : 2 ] != 'TL' ):
        f.seek(5 * 4 * recl); #reposition to end of record 5
        Pos.s.x     = unpack(str(Nsx)+'f',  f.read(Nsx*4))
        
        f.seek( 6 * 4 * recl); #reposition to end of record 6
        Pos.s.y     = unpack(str(Nsy) + 'f', f.read(Nsy*4))
    else:   # compressed format for TL from FIELD3D
        f.seek(5 * 4 * recl, -1 ); #reposition to end of record 5
        Pos.s.x     = f.read(2,    'float32' )
        Pos.s.x     = linspace( Pos.s.x( 1 ), Pos.s.x( end ), Nsx )
        
        f.seek(6 * 4 * recl, -1 ); #reposition to end of record 6
        Pos.s.y     = f.read(2,    'float32' )
        Pos.s.y     = linspace( Pos.s.y( 1 ), Pos.s.y( end ), Nsy )

    f.seek(7 * 4 * recl); #reposition to end of record 7
    Pos.s.depth = unpack(str(Nsd)+'f', f.read(Nsd*4))

    f.seek(8 * 4 * recl); #reposition to end of record 8
    Pos.r.depth = unpack(str(Nrd) + 'f', f.read(Nrd*4))

    f.seek(9 * 4 * recl); #reposition to end of record 9
    Pos.r.range = unpack(str(Nrr) + 'f',f.read(Nrr*4))
    # Pos.r.range = Pos.r.range';   # make it a row vector

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
        
        xdiff = abs( Pos.s.x - xs * 1000. )
        [ holder, idxX ] = min( xdiff )
        ydiff = abs( Pos.s.y - ys * 1000. )
        [ holder, idxY ] = min( ydiff )
        
        # show the source x, y that was found to be closest
        # [ Pos.s.x( idxX ) Pos.s.y( idxY ) ]
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
    return [ title, PlotType, freqVec, atten, Pos, pressure ] 


def read_shd ( *varargin ):

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
            [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] = read_shd_bin( filename )
        if len(varargin) ==  2:
            [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] = read_shd_bin( filename, freq )
        if len(varargin) ==  3:
            [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] = read_shd_bin( filename, xs, ys )
    elif FileType ==  'shdmat':  # Shade function mat file
        loadmat( filename )

        # has a specific source xs, ys been given?
        if not np.isnan( xs ):
            xdiff = abs( Pos.s.x - xs * 1000. )
            [ holder, idxX ] = min( xdiff )
            ydiff = abs( Pos.s.y - ys * 1000. )
            [ holder, idxY ] = min( ydiff )

            # extract the appropriate source index
            pressureT = pressure[ idxX, idxY, :, :, :, : ]
            pressure = np.reshape( pressureT, [ Pos.Nsd, Pos.Nrd, Pos.Ntheta, Pos.Nrr ] )

        # has a specific frequency been given?
        if not np.isnan( freq ):
            freqdiff = abs( freqVec - freq )
            [ holder, ifreq ] = min( freqdiff )
            # extract the appropriate source index
            pressure = pressure[ ifreq, 1, :, : ]

    elif FileType ==  'asc': # ascii format
        [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] = read_shd_asc( filename )

    elif FileType ==  'grnmat':   # Green's function mat file
        loadmat( filename )
        Pos.r.range = Pos.r.range.T;   # make it a column vector to match read_shd_bin

    elif FileType ==  'RAM':
        [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] = read_ram_tlgrid
    else:
        raise ValueError( 'Unrecognized file extension' )

    # clean up PlotTitle by taking only the part up inside the quotes
    # nchars = strfind( PlotTitle, '''' );   # find quotes
    # PlotTitle = [ PlotTitle( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) ]

    return [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] 
