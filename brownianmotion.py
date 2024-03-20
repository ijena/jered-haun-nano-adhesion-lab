import matplotlib.pyplot as plt
import math
import random
import numpy as np
from scipy.stats import norm

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
    velocity = 10  # initial velocity of the particle is assumed to be 0 nm/ns CHANGED
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
    average_particle_position = [
        x_position
    ]  # create a list to store the average particle position at every instant
    average_particle_velocity = [
        velocity
    ]  # create a list to store the average particle velocity at every instant
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
    last_position = 500
    last_velocity = 10  # CHANGED
    average_position = 0  # variable to store average position of the particle
    average_velocity = 0  # variable to store average velocity of the particle
    for time in range(1, total_time + 1):
        last_position = x_position
        last_velocity = velocity

        # calculating random sigma value for gaussian between -3sigma_position and 3sigma_position
        random_sigma_position = np.random.uniform(
            -3 * sigma_position, 3 * sigma_position
        )
        # equation to calculate new position using equation 5 in Hammer English paper
        # absolute value of random_sigma_position is taken because sigma(standard deviation) cannot be negative

        x_position = (
            x_position
            + (c1 * time_interval * velocity)
            + (c2 * pow(time_interval, 2) * K)
            + gaussian(last_position, np.abs(random_sigma_position))
        )
        # calculating random sigma value for gaussian between 0 and 3sigma_velocity

        random_sigma_velocity = np.random.uniform(0, 3 * sigma_velocity)
        # equation to calculate new velocity using equation 5 in Hammer English paper
        # absolute value of random_sigma_velocity is taken because sigma(standard deviation) cannot be negative
        velocity = (
            (c0 * velocity)
            + (c1 * time_interval * K)
            + gaussian(last_velocity, np.abs(random_sigma_velocity))
        )
        # making sure the particle stays within bounds
        x_position = x_position % (upper_bound - lower_bound)
        particle_positions.append(x_position)
        # calculate average position at that instant and add to the list to plot later
        average_position = np.mean(particle_positions)
        average_particle_position.append(average_position)
        # print(particle_positions)
        particle_velocity.append(velocity)
        particle_time.append(time)
        # calculate average velocity at that instant and add to the list to plot later
        average_velocity = np.mean(particle_velocity)
        average_particle_velocity.append(average_velocity)

    # print(particle_positions)
    # print(particle_velocity)
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    # TO DO - plot it around position = initial_position and velocity = initial_velocity
    # Plotting particle position vs time and the average position at each instant

    axs[0, 0].plot(particle_time, particle_positions, label="Position")
    axs[0, 0].plot(particle_time, average_particle_position, label="Average position")

    axs[0, 0].set_ylabel("Particle Position (nm)")
    axs[0, 0].set_xlabel("Time passed since start of simulation (ns)")
    axs[0, 0].legend()
    axs[0, 0].set_title(
        f"Position vs Time for Brownian Motion on a {line_length} nm Line"
    )

    # Plotting particle velocity vs time
    axs[0, 1].plot(particle_time, average_particle_velocity, label="Average velocity")
    axs[0, 1].plot(particle_time, particle_velocity, label="Velocity")
    axs[0, 1].set_ylabel("Velocity (nm/s)")
    axs[0, 1].set_xlabel("Time (ns)")
    axs[0, 1].legend()
    axs[0, 1].set_title(
        f"Velocity vs Time for Brownian Motion on a {line_length} nm Line"
    )

    # Calculate mean and standard deviations for plotting a normal distribution overlay
    mean_position = np.mean(particle_positions)
    position_std_dev = np.std(particle_positions)
    # Create a normal distribution and plot it
    position_normal_distribution = norm.pdf(
        np.linspace(0, 1000), mean_position, position_std_dev
    )

    # Plotting Histogram of particle_positions
    axs[1, 0].hist(particle_positions, bins=50, label="Positions")
    axs[1, 0].plot(
        np.linspace(0, 1000),
        position_normal_distribution,
        label="Normal distribution overlay",
    )
    axs[1, 0].set_ylabel("Frequency")
    axs[1, 0].set_xlabel("Position(nm)")
    axs[1, 0].set_title(
        f"Histogram for particle position in Brownian Motion on a {line_length} nm Line"
    )
    axs[1, 0].legend()

    # Plotting Histogram of particle_velocities
    axs[1, 1].hist(particle_velocity, label="Velocity", color="orange")
    axs[1, 1].set_ylabel("Frequency")
    axs[1, 1].set_xlabel("Velocity (nm/s)")
    axs[1, 1].set_title(
        f"Histogram for particle velocity in Brownian Motion on a {line_length} nm Line"
    )

    plt.show()


def gaussian(last_position, sigma):
    # generate a value in a gaussian distribution based on the mean being the previous position/velocity
    x = np.random.normal(last_position, sigma)
    return x


if __name__ == "__main__":
    brownian_motion_simulation()
