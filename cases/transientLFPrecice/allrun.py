
import os
from pyFSI.execution.caseManagement import *

# Define case name and location
caseName = "transientLFPrecice"


# Clean and run
cleanCase(caseName, solid=False)
case = runCase(caseName)

# Run the plots file
#os.system("python3 plots.py 1")
exec(open("plots.py").read()) # Equivalent



