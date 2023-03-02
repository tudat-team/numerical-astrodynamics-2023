###########################################################################
#
# # Numerical Astrodynamics 2022/2023
#
# # Assignment 1, Question 5 - Propagation Settings
#
###########################################################################


''' 
Copyright (c) 2010-2020, Delft University of Technology
All rights reserved

This file is part of Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

import os

import numpy as np
from matplotlib import pyplot as plt

from tudatpy.io import save2txt
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.astro import element_conversion


class ThrustGuidance:

    def __init__(self,
                 maximum_thrust: float, # Maximum thrust value that is to be used
                 true_anomaly_threshold: float, # Limiting value of true anomaly before and after node at which thrust should be on/off
                 bodies: environment.SystemOfBodies):
        self.maximum_thrust = maximum_thrust
        self.true_anomaly_threshold = true_anomaly_threshold
        self.bodies = bodies

    def compute_thrust_direction(self, current_time: float):

        # Check if computation is to be done. NOTE TO STUDENTS: ALL CALCULATION OF THRUST DIRECTION MUST BE INSIDE
        # THE FOLLOWING BLOCK
        if( current_time == current_time ):

            # Retrieve current JUICE Cartesian state w.r.t. Ganymede from environment
            current_cartesian_state = self.bodies.get_body( 'JUICE' ).state - self.bodies.get_body( 'Ganymede' ).state
            gravitational_parameter = self.bodies.get_body( 'Ganymede' ).gravitational_parameter
            current_keplerian_state = element_conversion.cartesian_to_keplerian( current_cartesian_state, gravitational_parameter )

            # Compute and return current thrust direction (3x1 vector)
            thrust_direction = XXXX
            XXXX

            # Here, the direction of the thrust (in a frame with inertial orientation; same as current_cartesian_state)
            # should be returned as a numpy unit vector (3x1)
            return XXXX

        # If no computation is to be done, return zeros
        else
            return np.zeros([3,1])
    def compute_thrust_magnitude(self, current_time: float):

        # Check if computation is to be done. NOTE TO STUDENTS: ALL CALCULATION OF THRUST MAGNITUDE MUST BE INSIDE
        # THE FOLLOWING BLOCK
        if( current_time == current_time ):

            # Retrieve current JUICE Cartesian  and Keplerian state w.r.t. Ganymede from environment
            current_cartesian_state = self.bodies.get_body( 'JUICE' ).state - self.bodies.get_body( 'Ganymede' ).state
            gravitational_parameter = self.bodies.get_body( 'Ganymede' ).gravitational_parameter
            current_keplerian_state = element_conversion.cartesian_to_keplerian( current_cartesian_state, gravitational_parameter )

            # Compute and return current thrust magnitude (scalar)
            thrust_magnitude = XXXX
            XXXX

            # Here, the value of the thrust magnitude (in Newtons, as a single floating point variable), should be returned
            return XXXX
        # If no computation is to be done, return zeros
        else:
            return 0.0




# Retrieve current directory
current_directory = os.getcwd()

# # student number: 1244779 --> 1244ABC
A = XXXX
B = XXXX
C = XXXX

simulation_start_epoch = 35.4 * constants.JULIAN_YEAR + A * 7.0 * constants.JULIAN_DAY + B * constants.JULIAN_DAY + C * constants.JULIAN_DAY / 24.0
simulation_end_epoch = simulation_start_epoch + 344.0 * constants.JULIAN_DAY / 24.0

###########################################################################
# CREATE ENVIRONMENT ######################################################
###########################################################################

# Load spice kernels.
spice.load_standard_kernels()
spice.load_kernel( current_directory + "/juice_mat_crema_5_1_150lb_v01.bsp" );

# Create settings for celestial bodies
bodies_to_create = XXXX
global_frame_origin = 'Ganymede'
global_frame_orientation = 'ECLIPJ2000'
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create, global_frame_origin, global_frame_orientation)

# Create environment
bodies = environment_setup.create_system_of_bodies(body_settings)

###########################################################################
# CREATE VEHICLE ##########################################################
###########################################################################

# Create vehicle object
bodies.create_empty_body( 'JUICE' )


###########################################################################
# CREATE THRUST MODEL #####################################################
###########################################################################

# Create thrust guidance object (e.g. object that calculates direction/magnitude of thrust)
thrust_magnitude = XXXX
true_anomaly_threshold = XXXX
thrust_guidance_object = ThrustGuidance( thrust_magnitude, true_anomaly_threshold, bodies )

# Create engine model (default JUICE-fixed pointing direction) with custom thrust magnitude calculation
constant_specific_impulse = XXXX
thrust_magnitude_settings = (
    propagation_setup.thrust.custom_thrust_magnitude_fixed_isp(
        thrust_guidance_object.compute_thrust_magnitude,
        constant_specific_impulse ) )
environment_setup.add_engine_model(
    'JUICE', 'MainEngine', thrust_magnitude_settings, bodies )

# Create vehicle rotation model such that thrust points in required direction in inertial frame
thrust_direction_function = thrust_guidance_object.compute_thrust_direction
rotation_model_settings = environment_setup.rotation_model.custom_inertial_direction_based(
    thrust_direction_function,
    "JUICE-fixed",
    "ECLIPJ2000" )
environment_setup.add_rotation_model( bodies, "JUICE", rotation_model_settings)


###########################################################################
# CREATE ACCELERATIONS ####################################################
###########################################################################

# Define bodies that are propagated, and their central bodies of propagation.
bodies_to_propagate = ['JUICE']
central_bodies = ['Ganymede']

# Define accelerations acting on vehicle.
acceleration_settings_on_vehicle = dict(
    XXXX
)

# Create global accelerations dictionary.
acceleration_settings = {'JUICE': acceleration_settings_on_vehicle}

# Create acceleration models.
acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, bodies_to_propagate, central_bodies)

###########################################################################
# CREATE PROPAGATION SETTINGS #############################################
###########################################################################

# Define initial state.
system_initial_state = spice.get_body_cartesian_state_at_epoch(
    target_body_name='JUICE',
    observer_body_name='Ganymede',
    reference_frame_name='ECLIPJ2000',
    aberration_corrections='NONE',
    ephemeris_time = simulation_start_epoch )

# Define required outputs
dependent_variables_to_save = XXXX

# Create numerical integrator settings.
fixed_step_size = 10.0
integrator_settings = propagation_setup.integrator.runge_kutta_4(
    fixed_step_size
)

# Create translational propagation settings.
termination_settings = propagation_setup.propagator.time_termination( simulation_end_epoch )
translational_propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    system_initial_state,
    simulation_start_epoch,
    integrator_settings,
    termination_settings,
    output_variables = dependent_variables_to_save
)

# Create mass propagator settings
XXXX

# Create combined mass and translational dynamics propagator settings
XXXX
propagator_settings = ...


###########################################################################
# PROPAGATE ORBIT #########################################################
###########################################################################

# Create simulation object and propagate dynamics.
dynamics_simulator = numerical_simulation.create_dynamics_simulator(
    bodies, propagator_settings )

# Retrieve all data produced by simulation
propagation_results = dynamics_simulator.propagation_results

# Extract numerical solution for states and dependent variables
state_history = propagation_results.state_history
dependent_variables = propagation_results.dependent_variable_history

###########################################################################
# SAVE RESULTS ############################################################
###########################################################################

save2txt(solution=state_history,
         filename='JUICEPropagationHistory_Q1.dat',
         directory='./'
         )

save2txt(solution=dependent_variables,
         filename='JUICEPropagationHistory_DependentVariables_Q1.dat',
         directory='./'
         )

###########################################################################
# PLOT RESULTS ############################################################
###########################################################################

# Extract time and Kepler elements from dependent variables
kepler_elements = np.vstack(list(dependent_variables.values()))
time = dependent_variables.keys()
time_days = [ t / constants.JULIAN_DAY - simulation_start_epoch / constants.JULIAN_DAY for t in time ]








