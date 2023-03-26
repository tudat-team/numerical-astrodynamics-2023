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
# RUN CODE FOR QUESTION 5 #################################################
###########################################################################

# Create body objects
bodies = create_simulation_bodies()

# Create Lambert arc state model
lambert_arc_ephemeris = get_lambert_problem_result(bodies, target_body, departure_epoch, arrival_epoch)


# Set arc length
number_of_arcs = 10
arc_length = XXXX

for arc_index in range(number_of_arcs)

    # Compute start and end time for current arc
    current_arc_initial_time = XXXX
    current_arc_final_time = XXXX

    # Get propagator settings for perturbed forward arc
    arc_initial_state = lambert_arc_ephemeris.cartesian_state(current_arc_initial_time)
    propagator_settings = get_perturbed_propagator_settings(bodies, arc_initial_state, current_arc_initial_time, current_arc_final_time)

    ###########################################################################
    # PROPAGATE NOMINAL TRAJECTORY AND VARIATIONAL EQUATIONS ##################
    ###########################################################################

    sensitivity_parameters = get_sensitivity_parameter_set(
        propagator_settings, bodies )
    variational_equations_simulator = numerical_simulation.create_variational_equations_solver(
        bodies, propagator_settings, sensitivity_parameters)

    state_transition_result = variational_equations_simulator.state_transition_matrix_history
    nominal_integration_result = variational_equations_simulator.state_history

    # Computer arc initial state before applying variations
    initial_epoch = list(state_transition_result.keys())[0]
    original_initial_state = nominal_integration_result[initial_epoch]

    ###########################################################################
    # START ANALYSIS ALGORITHM FOR QUESTION 4 #################################
    ###########################################################################

    # This vector will hold the maximum permitted initial state perturbations for which the linearization 
    # is valid (for the current arc. The vector is initialized to 0, and each of its 6 entries is computed 
    # in the 6 iterations of the coming for loop (that runs over the iteration variable 'entry')
    permitted_perturbations = np.array([0, 0, 0, 0, 0, 0])

    # Iterate over all initial state entries
    for entry in range(6):

        # Define (iterative) algorithm to compute current entry of 'permitted_perturbations'
        # General structure: define an initial state perturbation (perturbed_initial_state variable),
        # compute epsilon_x (see assignment), and iterate your algorithm until convergence.

        while XXXX:

            # Reset propagator settings with perturbed initial state
            perturbed_initial_state = XXXX
            propagator_settings.initial_states = perturbed_initial_state

            XXXX

            # Compute epsilon_x
            epsilon_x = XXXX

        permitted_perturbations[entry] = XXXX
