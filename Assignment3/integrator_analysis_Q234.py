'''
Copyright (c) 2010-2020, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

#  This file should run the code for questions 2, 3 and 4.
#
#  The simulations for questions 2 and 3 are largely similar, with the exceptions that:#
#   * Question 2 uses perturbations, question 3 does not
#   * Question 2 requires a numerical benchmark, question 3 does not
#
#  To run question 4, toggle the run_question_4 variable to True. Running code for this questions modifies
#  the initial time, and limits the analysis to the flyby phase. To run the code for question 4, write the code for the
#  get_closest_approach_time function in the integrator_analysis_helper_functions file first.
#
#  Note that for question 4, you should (can) use both the perturbed and unperturbed dynamics, which are automatically
#  performed if you finish the script and set run_question_4 to True

from integrator_analysis_helper_functions import *

spice_interface.load_standard_kernels()
bodies = create_bodies( )

run_question_4 = False
if( run_question_4 ):
    number_of_iterations = 2
else:
    number_of_iterations = 1

# Define list of integrator tolerances
integration_tolerances = [1.0E-12, 1.0E-10, 1.0E-8, 1.0E-6]

# Run code for Q2 and Q3 (run_flyby_from_closest_approach=0) and optionally for Q4 (run_flyby_from_closest_approach=1)
for run_flyby_from_closest_approach in range(number_of_iterations):

    # Only run flyby phase for Q4
    if run_flyby_from_closest_approach == 1:
        number_of_phases = 1
    else:
        number_of_phases = 2

    # Iterate over the mission phases
    for current_phase in range( number_of_phases ):

        # Create initial state and time
        current_phase_start_time = ...
        current_phase_end_time = ...

        # If running Q4, modify initial/final times
        if run_flyby_from_closest_approach == 1:
            ...

        # Define central body of propagation 
        current_central_body = central_bodies_per_phase[ current_phase ]

        # Define termination conditions (enforce exact termination time)
        termination_condition = propagation_setup.propagator.time_termination(
            current_phase_end_time,
            terminate_exactly_on_final_condition=True)

        # Create acceleration models for perturbed and unperturbed case
        perturbed_acceleration_models = get_perturbed_accelerations( current_central_body, bodies)
        unperturbed_acceleration_models = get_unperturbed_accelerations( current_central_body, bodies)

        # Create propagator settings for perturbed and unperturbed case
        perturbed_propagator_settings = ...
        unperturbed_propagator_settings = ...

        # Define integrator settings for benchmark
        benchmark_integrator_settings = ...

        # Propagate benchmark dynamics
        benchmark_dynamics_simulator = ...

        # Create interpolator for benchmark results
        interpolator_settings = interpolators.lagrange_interpolation( 8 )
        benchmark_interpolator = interpolators.create_one_dimensional_interpolator(
            benchmark_dynamics_simulator.state_history, interpolator_settings )
        
        # Perform integration of dynamics with different tolerances
        for current_tolerance in integration_tolerances:
            
            # Define integrator step settings
            initial_time_step = 10.0
            minimum_step_size = 1.0E-16
            maximum_step_size = np.inf
            
            # Retrieve coefficient set
            coefficient_set = ...
            
            # Create variable step-size integrator settings
            integrator_settings = ...
            
            # Define output file name
            file_output_identifier = "Iteration_" + str(run_flyby_from_closest_approach) + "tolerance-index" +\
                                     str(integration_tolerances.index(current_tolerance)) + \
                                     "_phase-index" + str(current_phase)
            
            # Propagate dynamics for perturbed and unperturbed case
            perturbed_dynamics_simulator = ...
            unperturbed_dynamics_simulator = ...

            write_propagation_results_and_benchmark_difference_to_file(
                    perturbed_dynamics_simulator,
                    file_output_identifier,
                    benchmark_interpolator)

            write_propagation_results_and_analytical_difference_to_file(
                    unperturbed_dynamics_simulator,
                    file_output_identifier + "_unperturbed",
                    bodies.get_body(current_central_body).gravitational_parameter)
    
