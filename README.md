# pyFSI v0.2.4

This is the second version of pyFSI, a code for solving fluid-structure interaction problems.

Currently we have three FSI models: Tosi's linear Leakage Flow, Saravia's nonlinear Leakege Flow and preciceCoupling.

## Change List
* New boundary conditions for simple supported beams. 
* Boundary conditions are now specified through tables.


## Dependencies
* Python 3
* ffmpeg for rendering video files `sudo apt get install ffmpeg
* preCICE

## Installation
* Clone the repository to your computer.
* Edit your .bashrc and change the file limit adding `ulimit -n 2048`. For parametric analysis, you may need to set this value to an even larger number; for example, for 1000 parametric analysis `ulimit -n 102400` is required.
## Usage

Create a folder in the cases directory containing the .json input file and an allrun.py which executes the solver. 

For shell interactive mode run `python3 -i allrun.py`. To return to the shell after running: For interactive mode run `python3 allrun.py`

Follow the examples in the case directory. Currently we support three type of analysis: eigenvalue extraction, nonlinear transient and preCICE-Calculix-pyFSI transient. 

## Errors

### Limit of open files exceeded. 
This is an error that is often generated during large parametric runs. Somehow even when pyFSI opens the files as a+, a new file is opened. 

#### Solution
Execute in a shell `ulimit -n 2048, or increase the limits until the errors disappear. 

Note that closing the file after creation is expensive, so it is preferred to let it opened. 

### Running from PyCharm allrun.py cannot find pyFSI
This does not happen if we run from the command line as `python3 allrun.py`

#### Solution
Start PyCharm from a terminal in which .bashrc has added the folder where pyFSI is located to PYTHOPATH. 


## Tickets
* Fix the file limit problem.  


## Installation of preCICE and Calculix
### Precice core library
This is the main C++ library. You can find the installation steps in https://www.precice.org/installation-overview.html

`wget https://github.com/precice/precice/releases/download/v2.2.0/libprecice2_2.2.0_focal.deb
sudo apt install ./libprecice2_2.2.0_focal.deb` 

### Python bindings
The Python bindings are needed by the pyFSI adapter (who couples the fields). The installation procedure is at https://www.precice.org/installation-bindings-python.
For short `pip3 install --user pyprecice` 

### Calculix
The full installation procedure can be found at https://www.precice.org/adapter-calculix-get-calculix.html.
A brief summary is given next.

#### Install Calculix dependecies
`sudo apt install libarpack2-dev libspooles-dev libyaml-cpp-dev`

#### Get Calculix 
`cd /opt
sudo wget http://www.dhondt.de/ccx_2.16.src.tar.bz2
sudo tar xvjf ccx_2.16.src.tar.bz2 
sudo rm ccx_2.16.src.tar.bz2`

#### Get the adapter and modify it according to your system
These are the commands to get the adapter and modify the Makefile to point to the right location of Calculix (/opt) instead of home. 
`cd /opt
sudo wget https://github.com/precice/calculix-adapter/archive/master.zip 
sudo unzip master.zip 
sudo rm master.zip
cd calculix-adapter-master
sudo gedit Makefile
`
Change the location of Calculix (variable CCX)

#### Compile Calculix with the adapter
`cd /opt/calculix-adapter-master
sudo make -j 4`

#### Check
You should now have a new executable ccx_preCICE in the bin/ folder of the adapter. You may move this file to a path known by your system, or add this to your PATH (careful when doing this!).

`cd /opt/calculix-adapter-master/bin
sudo cp ccx_preCICE /usr/bin/.`

Probably it would be neccessary to move calculix-adapter-master to the Calculix directory
and also rename it to calculix-adapter.

`sudo mv opt/calculix-adapter-master opt/calculix-adapter
sudo mv opt/calculix-adapter /opt/CalculiX/.`

### Run precice case
It is often necessary to run the case twice from the same terminal without closing the additional terminal opened by pyFSI. We are now sure which is the problem.


