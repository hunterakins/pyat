Test specific documentation for PyAT
Hunter Akins
2020


------ 
Acoustic Toolbox install 

To run the tests you need a working installation of the Acoustics Toolbox (AT).
Once you build it, you should be able to go to the folder where you downloaded AT and go to each
model subdirectory and there will be a file there called Model.exe
For example:
>> cd at_project_directory
>> cd Kraken
>> ls
Blah blah blah lots of files and one of them is called "kraken.exe" if you did everything right.

Similarly for Bellhop and Scooter.

Python needs to be able to find the binaries. 
A nice solution to this problem is to add the location of the binaries (the .exe files) to your PATH environment variable. This works on Linux, Mac, and Windows as far as I know. 

-------
Adding binaries to your PATH so that you can call them from Python conveniently

What I did was created a subfolder in the AT directory for the binaries 
Suppose at_project_directory is the directory where you downloaded and compiled AT
(for me it's /home/hunter/Downloads/at)

>> cd at_project_directory
>> mkdir bin
>> cp Kraken/kraken.exe bin
>> cp Kraken/krakenc.exe bin
>> cp Bellhop/bellhop.exe bin
>> cp Bellhop/bellhop3d.exe bin
>> cp KrakenField/field.exe  bin
>> cp Scooter/scooter.exe bin
>> cp Scooter/fields.exe bin

Now, for example in my case, if I type
>> ls /home/hunter/Downloads/at/bin
I see:
bellhop.exe  field.exe  fields.exe  krakenc.exe  kraken.exe  scooter.exe

So at_project_directory/bin contains all the binaries you need for pyat

Now simply add at_project_directory/bin to your environment path. To do so: google "how to add new locations to my PATH variable in [windows, mac, whatever]).
Once you get to the part where you add the variable, just add "at_project_directory/bin"
Use absolute paths (e.g. on Linux /home/your_name/blahblah)

On Linux:
Open up .zshrc or .bashrc
vim ~/.zshrc

Add this line
export PATH="at_project_directory/bin:$PATH"
Make sure at_project_directory is your absolute path to the AT code

Restart your shell
source ~/.zshrc 
or source ~/.bashrc

Now if you type in field.exe to your shell, it will know where to look to run it

---------
Running the AT programs within Python

I use os.system to pass system calls as strings
There are other modules. I think subprocess might be good if you want to do threading.
os.system is sufficient for me

Example
>> from os import system
>> system('field.exe file')



----------
Basic structure of the tests

The tests mostly consist of a lot of lines of code where I configure the environment. 
I have to set the SSP, the source depth, number of sources, receiver positions, bottom properties.
Once these variables are set, you pass them to write_env
write_env will create a .env file in the test directory
Then you run a model
system('model_name env_file_name')

I think 2D bellhop is the best commented. It's 2D because it's a range dependent SSP


