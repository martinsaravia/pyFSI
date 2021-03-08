# Run the case
import os
from pyFSI.execution.caseManagement import *

# Define case name and location
caseName = "eigenLFBeam"


# Clean and run
cleanCase(caseName)
case = runCase(caseName)

# Run the plots file
#os.system("python3 plots.py 1")
exec(open("plots.py").read()) # Equivalent
