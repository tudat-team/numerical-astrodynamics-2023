'''
Copyright (c) 2010-2020, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''


from integrator_analysis_helper_functions import *

# Load spice kernels.
spice_interface.load_standard_kernels()

# Create the bodies for the numerical simulation
bodies = create_bodies( )

# Define list of step size for integrator to take
step_sizes = ...

# Iterate over phases
for current_phase in range( len(central_bodies_per_phase )):
    
    # Create initial state and time
    current_phase_start_time = initial_times_per_phase[ current_phase ]
    current_phase_end_time = current_phase_start_time + propagation_times_per_phase[ current_phase ]

    # Define current centra
    current_central_body = central_bodies_per_phase[ current_phase ]

    # Retrieve JUICE initial state
    initial_state = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name="JUICE",
        observer_body_name=current_central_body,
        reference_frame_name=global_frame_orientation,
        aberration_corrections="NONE",
        ephemeris_time=current_phase_start_time
    )
    
    # Retrieve acceleration settings without perturbations
    acceleration_models = get_unperturbed_accelerations( current_central_body, bodies)

    # Define propagator settings
    propagator_settings = ...
    
    # Iterate over step size
    for step_size in step_sizes:
        
        # Define integrator settings
        integrator_settings = get_fixed_step_size_integrator_settings(current_phase_start_time, step_size)

        # Propagate dynamics
        dynamics_simulator = ...
        state_history = dynamics_simulator.state_history
        
        # Compute difference w.r.t. analytical solution to file
        central_body_gravitational_parameter = bodies.get_body( current_central_body ).gravitational_parameter
        keplerian_solution_difference = get_difference_wrt_kepler_orbit( 
            state_history, central_body_gravitational_parameter)
        
        # Write results to files
        file_output_identifier = "Q1_step_size" + str(step_size) + "_phase_index" + str(current_phase)   
        write_propagation_results_and_analytical_difference_to_file( 
            dynamics_simulator, file_output_identifier, bodies.get_body( current_central_body ).gravitational_parameter)
        
