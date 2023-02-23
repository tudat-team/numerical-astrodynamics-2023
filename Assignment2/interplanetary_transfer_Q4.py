''' 
Copyright (c) 2010-2020, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

from interplanetary_transfer_helper_functions import *

# Load spice kernels.
spice_interface.load_standard_kernels( )

# Define directory where simulation output will be written
output_directory = "./SimulationOutput/"

###########################################################################
# RUN CODE FOR QUESTION 4 #################################################
###########################################################################

rsw_acceleration_magnitude = [0, 0, 0]

# Create body objects
bodies = create_simulation_bodies()

# Create Lambert arc state model
lambert_arc_ephemeris = get_lambert_problem_result(bodies, target_body, departure_epoch, arrival_epoch)

###########################################################################
# RUN CODE FOR QUESTION 4b ################################################
###########################################################################

# Set start and end times of full trajectory
departure_epoch_with_buffer = XXXX
arrival_epoch_with_buffer = XXXX

# Solve for state transition matrix on current arc
variational_equations_solver = propagate_variational_equations(
    departure_epoch_with_buffer,
    arrival_epoch_with_buffer,
    bodies,
    lambert_arc_ephemeris,
    use_rsw_acceleration = True)

sensitivity_matrix_history = variational_equations_solver.sensitivity_matrix_history
state_history = variational_equations_solver.state_history
lambert_history = get_lambert_arc_history(lambert_arc_ephemeris, state_history)

# Compute low-thrust RSW acceleration to meet required final position
rsw_acceleration_magnitude = XXXX

# Propagate dynamics with RSW acceleration. NOTE: use the rsw_acceleration_magnitude as
# input to the propagate_trajectory function
dynamics_simulator = XXXX

###########################################################################
# RUN CODE FOR QUESTION 4e ################################################
###########################################################################

