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
# RUN CODE FOR QUESTION 3 #################################################
###########################################################################

# Create body objects
bodies = create_simulation_bodies()

# Create Lambert arc state model
lambert_arc_ephemeris = get_lambert_problem_result(bodies, target_body, departure_epoch, arrival_epoch)

##############################################################
# Compute number of arcs and arc length
number_of_arcs = 10
arc_length = XXXX

##############################################################

# Compute relevant parameters (dynamics, state transition matrix, Delta V) for each arc
for arc_index in range(number_of_arcs):

    # Compute initial and final time for arc
    current_arc_initial_time = XXXX
    current_arc_final_time = XXXX

    ###########################################################################
    # RUN CODE FOR QUESTION 3a ################################################
    ###########################################################################

    # Propagate dynamics on current arc (use propagate_trajecory function)
    XXXX

    ###########################################################################
    # RUN CODE FOR QUESTION 3c/d/e ############################################
    ###########################################################################
    # Note: for question 3e, part of the code below will be put into a loop
    # for the requested iterations

    # Solve for state transition matrix on current arc
    variational_equations_solver = propagate_variational_equations(current_arc_initial_time,
                                                                   current_arc_final_time, bodies,
                                                                   lambert_arc_ephemeris)
    state_transition_matrix_history = variational_equations_solver.state_transition_matrix_history
    state_history = variational_equations_solver.state_history
    lambert_history = get_lambert_arc_history(lambert_arc_ephemeris, state_history)

    # Get final state transition matrix (and its inverse)
    final_epoch = list(state_transition_matrix_history.keys())[-1]
    final_state_transition_matrix = state_transition_matrix_history[final_epoch]

    # Retrieve final state deviation
    final_state_deviation = state_history[final_epoch] - lambert_history[final_epoch]

    # Compute required velocity change at beginning of arc to meet required final state
    initial_state_correction = XXXX

    # Propagate with correction to initial state (use propagate_trajecory function),
    # and its optional initial_state_correction input
    dynamics_simulator = XXXX

