import matplotlib.pyplot as plt
import math
import random
import numpy as np

# for plots, plot position vs time and velocity vs time, around x = 0 or v = 0 ( 0 as in initial position which is middle of the box)
# histograms of velocity and position which should look like a normal distribution
# calculate sigma using equation 6 from hammer english
# use nanoseconds and nanometre

# TO DO add error handling for inputs
# def value_error_handling(user_input, data_type):
#    try:


def brownian_motion_simulation():

    line_length = 1000  # length of box in nanometre

    x_position = 500  # initial position of the particle at midpoint of the line which is at the 500 nm
    velocity = 0  # initial velocity of the particle is assumed to be 0 nm/ns
    total_time = 1000  # time of simulation in nanoseconds
    particle_radius = 100  # particle_radius in nanometre
    particle_mass = (1.05 / pow(10, 12)) * (
        (4.0 / 3.0) * math.pi * pow(particle_radius, 3)
    )  # mass = density * volume of polysterene. density is in nanogram/nm^3 and volume is in nm^3
    # lower_bound and upper_bound refer to the endpoints of the box
    lower_bound = 0
    upper_bound = 1000
    particle_positions = [
        x_position
    ]  # particle_positions is created to store all the positions of the particle to plot later. x_position is the initial position of the particle at time = 0 ns
    particle_velocity = [
        velocity
    ]  # particle_velocity is created to store all the velocities of the particle to plot later. initial velocity is assumed to be 0 nm/ns
    particle_time = [
        0
    ]  # particle_time is created to store the time passed since start of simulation. 0ns is the start of the simulation
    # check unit conversion
    viscosity_liquid = 1 / pow(10, 18)  # viscosity of water in ng/(nm.ns)
    inverse_viscous_relaxation_time = (
        3 * math.pi * viscosity_liquid * 2 * particle_radius
    ) / particle_mass  # calculating inverse viscous relaxation time using equation 4 from Hammer- English paper
    boltzmann_constant = 1.380649 * pow(
        10, -11
    )  # boltzmann constant in (nm^2 * ng)/(ns^2 * K)
    temperature = 310.2  # average body temperature in Kelvin
    time_interval = 1  # time_interval is taken as 1 ns
    # constants c0,c1 and c2 from equation 5 in Hammer English paper
    c0 = math.exp(-inverse_viscous_relaxation_time * time_interval)
    c1 = (1 - c0) / (inverse_viscous_relaxation_time * time_interval)
    c2 = (1 - c1) / (inverse_viscous_relaxation_time * time_interval)
    K = 0  # K in equations is assumed to be 0 currently
    # equation to calculate sigma_position using equation 6 in Hammer English paper
    sigma_position = math.sqrt(
        pow(time_interval, 2)
        * ((boltzmann_constant * temperature) / particle_mass)
        * (
            2
            - (
                1
                / (inverse_viscous_relaxation_time * time_interval)
                * (
                    3
                    - 4 * math.exp(-inverse_viscous_relaxation_time * time_interval)
                    + math.exp(-2 * inverse_viscous_relaxation_time * time_interval)
                )
            )
        )
    )
    # equation to calculate sigma_velocity using velocity 6 in Hammer English paper
    sigma_velocity = math.sqrt(
        ((inverse_viscous_relaxation_time * temperature) / particle_mass)
        * (1 - math.exp(-2 * inverse_viscous_relaxation_time * time_interval))
    )
    for time in range(1, total_time + 1):

        # equation to calculate new position using equation 5 in Hammer English paper
        x_position = (
            x_position
            + (c1 * time_interval * velocity)
            + (c2 * pow(time_interval, 2) * K)
            + gaussian(sigma_position)
        )
        # equation to calculate new velocity using equation 5 in Hammer English paper
        velocity = (c0 * velocity) + (c1 * time_interval * K) + gaussian(sigma_velocity)
        # the position of the particle is checked in a loop if it is within bounds
        while x_position < lower_bound or x_position > upper_bound:
            # if the position of particle goes out of the lower bound it continues from the upper bound as if the box has infinite length
            if x_position < lower_bound:
                out_of_bounds_position = lower_bound - x_position
                x_position = upper_bound - out_of_bounds_position
            # test
            # if the position of the particle goes out of the upper bound it continues from the lower bound as if the box has infinite length
            elif x_position > upper_bound:
                out_of_bounds_position = x_position - upper_bound
                x_position = lower_bound + out_of_bounds_position
        particle_positions.append(x_position)
        particle_velocity.append(velocity)
        particle_time.append(time)

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    # TO DO - plot it around position = initial_position and velocity = initial_velocity
    # Plotting particle position vs time
    axs[0, 0].plot(particle_time, particle_positions, label="Position")
    axs[0, 0].set_ylabel("Particle Position (nm)")
    axs[0, 0].set_xlabel("Time passed since start of simulation (ns)")
    axs[0, 0].set_title(
        f"Position vs Time for Brownian Motion on a {line_length} nm Line"
    )

    # Plotting particle velocity vs time
    axs[0, 1].plot(particle_time, particle_velocity, label="Velocity", color="orange")
    axs[0, 1].set_ylabel("Velocity (nm/s)")
    axs[0, 1].set_xlabel("Time (ns)")
    axs[0, 1].set_title(
        f"Velocity vs Time for Brownian Motion on a {line_length} nm Line"
    )

    # TO DO - add normal distribution curve to the histograms
    # Plotting Histogram of particle_positions
    axs[1, 0].hist(particle_positions, label="Positions")
    axs[1, 0].set_ylabel("Frequency")
    axs[1, 0].set_xlabel("Position(nm)")
    axs[1, 0].set_title(
        f"Histogram for particle position in Brownian Motion on a {line_length} nm Line"
    )

    # Plotting Histogram of particle_velocities
    axs[1, 1].hist(particle_velocity, label="Velocity", color="orange")
    axs[1, 1].set_ylabel("Frequency")
    axs[1, 1].set_xlabel("Velocity (nm/s)")
    axs[1, 1].set_title(
        f"Histogram for particle velocity in Brownian Motion on a {line_length} nm Line"
    )

    plt.show()


def gaussian(sigma):
    x = np.random.normal(500, 3 * sigma)  # bounds of x
    # use gaussian distribution to return random vector
    return (1 / (sigma * math.sqrt(2 * math.pi))) * math.exp(
        -0.5 * pow(((x - 0) / sigma), 2)
    )


if __name__ == "__main__":
    brownian_motion_simulation()
