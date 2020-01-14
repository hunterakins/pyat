The code is all in the pyat subdirectory

To install:

cd to a directory in your Python path

git clone https://github.com/hunter/pyat

That will install the code

To check it worked
Start python interpreter

>> python
>>>> import pyat
>>>> from pyat.pyat.readwrite import read_shd

Those are just some sample imports


------------------------
The code relies on a successful install of the acoustics toolbox
(https://oalib-acoustics.org/AcousticsToolbox/index_at.html)
There are instructions there how to build it. 
It usually requires a small tweak of the Makefile because of the fortran libraries on your local machine.
It's worth getting it running those it's very powerful. 

Once you build the acoustics toolbox, you get access to a set of binaries (with a dummy extension .exe even on Linux) that take in text configuration files.
In principal, you could never use matlab or python, and manually create environment files, field files, etc. and manually run the models on them e.g. >> /home/yourname/directorywherethebinariesare/kraken.exe manually_created_env_file.env
In practice though you want to plot and process the simulation outputs. 
Currently the at contains a matlab interface to reading and writing the env files and output files (described briefly below)
I basically just translated the matlab to Python so you can do the same things with Python (read the output to numpy arrays, etc.)

My pyat code is mostly designed around the reading and writing of the various input output files.
A brief summary follows. The more comprehensive documentation is found in the acoustic toolbox codebase
- .env files
    Contain information about source location, receiver positions, sound speed profile, bathymetry, and boundary types.
- .mod files
    Output of a Kraken run. Contain modal wavenumbers, mode shapes at source and rcvr positions, other info.
- .flp files
    Configuration for a "field" run. field is a program that will take in a set of grid positions and compute the green's function on the grid. 
    With Kraken, you run field after you run kraken since it relies on the .mod file to compute the field
    With Scooter, you run field after you run scooter since it relies on the .grn file to copmute the field
- .grn file
    output of Scooter run, contains the wavenumber domain depth green's function. 
- .shd files
    Output of a field run. Field will take a .mod file or a .grn file and produce the green's function at the positions specified
    Use read_shd to read in the shd file to a python numpy array (see tests)


--------------------------
What is where and how do I find it?
pyat folder contains code
readwrite is the translation of read_env, write_env, read_shd, read_modes, blah blah blah
env defines a number of objects used to represent the elements of the acoustic model
There are objects for Source, Dom (domain), Pos (source and domain together), Modes, and more
Pretty straightforward stuff.
I will be adding more documentation in the files in the next couple months. 





