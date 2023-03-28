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

spice_interface.load_standard_kernels()
bodies = create_bodies( )

# Define settings for orbit only, unperturbed only
current_phase = 1
central_body = "Ganymede"
bodies_to_integrate = ["JUICE"]
central_bodies = [ central_body ]
current_phase_start_time = initial_times_per_phase[current_phase]
acceleration_models = get_unperturbed_accelerations( central_body, bodies)

state_differences_rkf = np.zeros((6,20))
state_differences_euler = np.zeros((6,20))

step_per_run = 600.0

# Perform 20 individual steps
for i in range(20):
    
    # Compute initial time of current step
    current_start_time = ...
    
    # Compute initial state of current step
    initial_state = ...
    
    '''
    RKF
    '''

    # Define propagator settings, terminate after 300 s.
    time_step = 300.0
    termination_time = current_start_time + time_step    
    termination_settings = propagation_setup.propagator.time_termination(
        termination_time, terminate_exactly_on_final_condition=True )
    propagator_settings = propagation_setup.propagator.translational(
            central_bodies,
            acceleration_models,
            bodies_to_integrate,
            initial_state,
            termination_settings )
    
    # Get fixed step RKF78 integrator settings
    integrator_settings = ...
    
    # Propagate Dynamics
    dynamics_simulator = ...
    local_truncation_error = ...
    
    '''
    Euler
    '''
    


