# pyFSI v0.2.2

This is the second version of pyFSI, a code for solving fluid-structure interaction problems.

Currently we have three FSI models: Tosi's linear Leakage Flow, Saravia's nonlinear Leakege Flow and preciceCoupling.

## Installation
* Clone the repository to your computer. 

* Edit your .bashrc and change the file limit adding `ulimit -n 2048`

## Usage

Create a folder in the cases directory containing the .json input file and an allrun.py which executes the solver. 

For shell interactive mode run `python3 -i allrun.py`. To return to the shell after running: For interactive mode run `python3 allrun.py`

Follow the examples in the case directory. Currently we support three type of analysis: eigenvalue extraction, nonlinear transient and preCICE-Calculix-pyFSI transient. 

## Errors

### Limit of open files exceeded. 
This is an error that is often generated during large parametric runs. Somehow even when pyFSI opens the files as a+, a new file is opened. 

#### Solution
Execute in a shell `ulimit -n 2048, or increase the limits until the errors dissapear. 

Note that closing the file after creation is expensive, so it is preferred to let it opened. 


## Tickets
* Fix the file limit problem. 