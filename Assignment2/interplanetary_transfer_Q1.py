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
# RUN CODE FOR QUESTION 1 #################################################
###########################################################################

# Create body objects
bodies = create_simulation_bodies( )

# Create Lambert arc state model
lambert_arc_ephemeris = get_lambert_problem_result(bodies, target_body, departure_epoch, arrival_epoch)

# Create propagation settings and propagate dynamics
dynamics_simulator = propagate_trajectory( departure_epoch, arrival_epoch, bodies, lambert_arc_ephemeris,
                     use_perturbations = False)

# Write results to file
write_propagation_results_to_file(
    dynamics_simulator, lambert_arc_ephemeris, "Q1",output_directory)

# Extract state history from dynamics simulator
state_history = dynamics_simulator.state_history

# Evaluate the Lambert arc model at each of the epochs in the state_history
lambert_history = get_lambert_arc_history( lambert_arc_ephemeris, state_history )

