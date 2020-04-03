The Python code is all in the pyat subdirectory

To install:

cd to a directory in your Python path

git clone https://github.com/hunterakins/pyat

That will install the code

The dependencies are all fairly standard and should be included in an anaconda3 installation. 

You need git to download it this way. You can also go to the link above and click Download
and use your browser interface to download it. Just make sure the code is in your PYTHONPATH 
so you can import it. 

To check it worked
Start python interpreter

>> python
>>>> import pyat
>>>> from pyat.pyat.readwrite import read_shd

Those are just some sample imports to make sure the code is in a place that Python looks. 
If there are errors with this, you may try 

>>>> from pyat.readwrite import read_shd

If there are more errors, try
>>>> from os import sys
>>>> sys.path
[ '', '/home/your_name/research/code', '/home/your_name/research/models', '/home/your_name/anaconda3', '/home/your_name/pylib', '/home/your_name', '/home/your_name/anaconda3/lib/python37.zip', '/home/your_name/anaconda3/lib/python3.7', '/home/your_name/anaconda3/lib/python3.7/lib-dynload', '/home/your_name/anaconda3/lib/python3.7/site-packages'

Check that the project lives in one of the folders that is in the example list I provided above

------------------------
Acoustics Toolbox:

The code relies on a successful install of the Acoustics Toolbox
(https://oalib-acoustics.org/AcousticsToolbox/index_at.html)
There are instructions there how to build it. 
If it fails to build, it's likely because you need to install some specific fortran libraries on your computer. 
Just look carefully at the makefile output and find out what is breaking the build.
Then google how to install those fortran libraries on your computer. 
It's worth putting in effort to compile it because it's so useful.

Brief overview of the AT (for completeness):
Once you build the acoustics toolbox, you get access to a set of binaries (with a dummy extension .exe even on Linux) that take in text configuration files.
In principal, you could never use matlab or python, and manually create environment files, field files, etc. and manually run the models on them e.g. >> /home/yourname/directorywherethebinariesare/kraken.exe manually_created_env_file.env
This command will fied the text file *.env to the kraken program. 
If the .env file is properly created, Kraken will output a .mod file
In practice though you want to plot and process the simulation outputs. 
Since .mod is a binary file type, you need some programs to read it into a useful format.
Currently the AT contains a matlab interface for doing just that: reading and writing the env files and output files (overview of filetypes given below)
I basically just translated the matlab to Python so you can do the same things with Python (read the output to numpy arrays, etc., write an array of SSP values to an .env file and run the binaries)

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
- .arr files
    Arrival structure output from a BELLHOP run. 

Look through the documentation inside the AT for more information. 

--------------------------------
PyAT

What is where and how do I find it?

- pyat folder contains code (if you want to read or tweak)
    - readwrite.py is the translation of read_env, write_env, read_shd, read_modes, blah blah blah
    - env.py defines a number of objects used to represent the elements of the acoustic model
    There are objects for Source, Dom (domain), Pos (source and domain together), Modes, and more
    Pretty straightforward stuff.
- test folder contains example tests where I create environments in Python and run the models
    There is some more specific documentation there
- at_copy
    A copy of the acoustic toolbox that I use is there. That one is guaranteed to be compatible with my code. 

- docs
    You're already here ;)





