from pyFSI.execution.case import *
from pyFSI.post.prints import *
import sys

# Define case name
caseName = "eigenLFBeam"

# Create a log file
sys.stdout = io.Logger()

# Create the case object
case = MFSICase(caseName)

# Solve the case
case.solve()

# Run the plots file
exec(open("plots.py").read())

# Print some data
printArray(case.MFSI.K, "FSI Stiffness Matrix")
printNumbers(case.MFSI)


