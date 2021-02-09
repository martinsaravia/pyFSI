# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 07/12/2020
#    F    #     author: Martin Saravia
#    S    #     description: BE Beam class
#    I    #     return: solidModel object
# --------------------------------------------------------------------------- #
# Notes:
#   Staggered solver for fsi using preCICE and Calculix
#   a) The solid is solved using Calculix
# --------------------------------------------------------------------------- #

import importlib, shutil, os, pathlib
from solvers.solverBase import solverBase
import scipy.integrate as si
import numpy as np
import precice as prc


# Solver for transient FSI simulations
class preCICE(solverBase):
    def __init__(self, fsi):
        super().__init__(fsi)
        self.cleanCase()  # Clean the case for optimal linking

    def solve(self):
        interface = self.fsi.interface
        meshes = self.fsi.meshes
        flow = self.fsi.flow()

        # preCICE defines timestep size of solver via precice-config.xml
        precice_dt = interface.initialize()

        # Read initial data
        if interface.is_read_data_available():
            print("Reading available initial data ...")
            for mesh in meshes:
                mesh.readFields(interface)
            print("Available initial data read")

        print("--> Setting the flow intial conditions...")
        totalTime = precice_dt  # Total simluation time
        flow.setInitialConditions()
        state = flow.Q0

        # Write initial force data data
        if interface.is_action_required(prc.action_write_initial_data()):
            print("--> Writing available initial data ...")
            for mesh in meshes:
                mesh.writeFields(interface, flow)
            interface.mark_action_fulfilled(prc.action_write_initial_data())

        # Start the solution loop
        print("--> Starting the preCICE loop ...")
        while interface.is_coupling_ongoing():
            # When an implicit coupling scheme is used, checkpointing is required
            if interface.is_action_required(prc.action_write_iteration_checkpoint()):
                interface.mark_action_fulfilled(prc.action_write_iteration_checkpoint())

            # Solve the flow
            print("--> Solving the flow for time: ", totalTime)
            solution = si.solve_ivp(flow.rhs,
                                    [totalTime, totalTime+precice_dt],
                                    state,
                                    method='Radau')
            # Update the flow state
            print("--> Updating the flow object...")
            flow.update(totalTime, solution.y[:, -1])

            # Update the pressure distribution
            print("--> Updating the flow pressure distribution...")
            flow.updateForces(totalTime)

            # Write the force in the precice interface
            print("--> Writing the fluid force to the partner mesh...")
            for mesh in meshes:
                mesh.writeFields(interface, flow)

            # Advance in time
            print("--> Advancing the precice interface...")
            precice_dt = interface.advance(precice_dt)  # Execute the mapping

            # Read the displacement and velocities and write the the pyFSI boundary objects
            print("--> Reading the solid fields...")
            for mesh in meshes:
                mesh.readFields(interface)

            # Update the region geometry after updating the boundaries
            print("--> Updating the region dynamics...")
            for i, region in enumerate(self.fsi.flow().regions):
                region.update()

            # Check convergence and confirm advancing
            print("--> Cheking convergence...")
            if interface.is_action_required(prc.action_read_iteration_checkpoint()):  # i.e. not yet converged
                interface.mark_action_fulfilled(prc.action_read_iteration_checkpoint())
            else:  # converged, timestep complete
                self.fsi.flow().write()
                totalTime += precice_dt

        print("Exiting precice...")

        interface.finalize()

        return flow

    # Clean the case (neccessary for optimal linking betwen the participants)
    def cleanCase(self):
        print("--> Allcleaning...")
        for p in pathlib.Path(".").glob("P*.log"):
            p.unlink()
        try:
            shutil.rmtree("precice-run")
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))
