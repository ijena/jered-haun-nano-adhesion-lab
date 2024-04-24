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
    last_velocity = 0  # CHANGED
    average_position = 0  # variable to store average position of the particle
    average_velocity = 0  # variable to store average velocity of the particle
    for time in range(1, total_time + 1):
        last_position = x_position
        last_velocity = velocity

        # calculating random sigma value for gaussian between -3sigma_position and 3sigma_position
        random_sigma_position = np.random.normal(500, sigma_position)
        # equation to calculate new position using equation 5 in Hammer English paper
        # absolute value of random_sigma_position is taken because sigma(standard deviation) cannot be negative

        x_position = x_position + (
            (c1 * time_interval * velocity)
            + (c2 * pow(time_interval, 2) * K)
            + gaussian(last_position, np.abs(random_sigma_position))
        )
        # calculating random sigma value for gaussian between 0 and 3sigma_velocity

        random_sigma_velocity = np.random.normal(10**5, sigma_velocity)
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

    plt.plot(particle_time, particle_positions, label="Position")
    plt.plot(particle_time, average_particle_position, label="Average position")

    plt.ylabel("Particle Position (nm)")
    plt.xlabel("Time passed since start of simulation (ns)")
    plt.legend()
    plt.title(f"Position vs Time for Brownian Motion on a {line_length} nm Line")
    plt.show()
    # finding the sum of absolute positions at each second to analyze the position-time graph
    particle_position_analysis = [abs(particle_positions[0])]
    # calculate the cumulative sum of absolute positions and store it in particle_position_analysis
    for time in range(1, total_time + 1):
        particle_position_analysis.append(
            particle_position_analysis[time - 1] + abs(particle_positions[time])
        )
    # plotting cumulative sum of absolute positions with respect to time
    plt.plot(
        particle_time,
        particle_position_analysis,
        label="Cumulative sum of absolute positions",
    )
    plt.ylabel("Cumulative position (nm)")
    plt.xlabel("Time(ns)")
    plt.legend()
    plt.title("Cumulative sum of absolute positions at each instance of time")
    plt.show()

    # calculate slope of the cumulative sum graph at each instant
    cumulative_sum_position_slope = [0]
    for time in range(1, total_time + 1):
        cumulative_sum_position_slope.append(
            calculate_slope(
                time - 1,
                particle_position_analysis[time - 1],
                time,
                particle_position_analysis[time],
            )
        )

    # print(cumulative_sum_position_slope)
    # now we have to check if the slope at each instant = ndt where n is the number of dimensions, d is the diffusion coefficient and t is the time
    # Using Stokes-Einstein-Sutherland equation to find the diffusion coefficient https://en.wikipedia.org/wiki/Einstein_relation_(kinetic_theory)
    diffusion_coefficient = (boltzmann_constant * temperature) / (
        6 * math.pi * viscosity_liquid * particle_radius
    )
    position_slope_analysis = [0]
    for time in range(1, total_time + 1):
        position_slope_analysis.append(1 * diffusion_coefficient * time)
    print(1 * diffusion_coefficient * temperature)
    # print("actual", cumulative_sum_position_slope)
    # print("ideal", position_slope_analysis)
    # print(slope_analysis)
    # Plotting what the slope is and what it should be at every time instant for cumulative sum of particle position
    plt.plot(
        particle_time, cumulative_sum_position_slope, label="Actual slope of the graph"
    )
    plt.plot(particle_time, position_slope_analysis, label="Ideal slope of the graph")
    plt.ylabel("Slope(nm/ns)")
    plt.xlabel("Time (ns)")
    plt.title(
        "Ideal vs Real slope of the graph of cumulative sum of absolute positions of the particle at every instant"
    )
    plt.legend()
    plt.show()
    # Plotting particle velocity vs time
    plt.plot(particle_time, average_particle_velocity, label="Average velocity")
    plt.plot(particle_time, particle_velocity, label="Velocity")
    plt.ylabel("Velocity (nm/s)")
    plt.xlabel("Time (ns)")
    plt.legend()
    plt.title(f"Velocity vs Time for Brownian Motion on a {line_length} nm Line")
    plt.show()

    velocity_analysis = (
        []
    )  # list to store all the cumulative particle velocities for mean velocities from 10^1 to 10^10
    cumulative_sum_velocity_slope = (
        []
    )  # list to store the slope for all the cumulative sum of absolute velocity graphs
    for index in range(1, 11):
        velocity = 0  # set initial velocity to 0
        current_cumulative_velocity = [0]
        current_cumulative_velocity_slope = [
            0
        ]  # list to store the slope at every instant for a particular mean in the gaussian distribution
        # list to append cumulative velocity for that  particular mean
        for time in range(1, total_time + 1):
            # calculate velocity at every instant for the different means
            last_velocity = velocity

            random_sigma_velocity = np.random.normal(10**index, sigma_velocity)
            # equation to calculate new velocity using equation 5 in Hammer English paper
            # absolute value of random_sigma_velocity is taken because sigma(standard deviation) cannot be negative
            velocity = (
                (c0 * velocity)
                + (c1 * time_interval * K)
                + gaussian(last_velocity, np.abs(random_sigma_velocity))
            )
            # calculate the slope at that instant of time and add it to the list current_cumulative_velocity_slope
            current_cumulative_velocity_slope.append(
                calculate_slope(
                    time - 1, current_cumulative_velocity[-1], time, abs(velocity)
                )
            )
            # add the cumulative  absolute velocity of the particle at every instant in current_cumulative_velocity
            # add the previous cumulative velocity and the current cumulative velocity to find the new cumulative velocity

            current_cumulative_velocity.append(
                current_cumulative_velocity[-1] + abs(velocity)
            )
        # save the current_cumulative_velocity for that particular mean velocity in velocity_analysis
        velocity_analysis.append(current_cumulative_velocity)
        # print(current_cumulative_velocity_slope)
        cumulative_sum_velocity_slope.append(current_cumulative_velocity_slope)

    # plotting graphs of the velocity analysis
    # ERROR Watch, anything over this range cannot be plotted, is it because of infinite values?
    for index in range(0, 5):
        plt.plot(
            particle_time,
            velocity_analysis[index],
            label=f"Mean = 10 ^{index+1} nm/ns ",
        )
    plt.ylabel("Velocity(nm/ns)")
    plt.xlabel("Time (ns)")
    plt.title(
        "Cumulative velocity of particles with different means of the gaussian distribution "
    )
    plt.legend()
    plt.show()

    # plot the actual and ideal slopes for velocity analysis
    # TO DO plot ideal slopes for vector analysis
    for index in range(0, 10):
        plt.plot(
            particle_time,
            cumulative_sum_velocity_slope[index],
            label=f"Actual slope for mean = 10^{index+1} nm/ns",
        )
    plt.ylabel("Values of the slope (nm/ns)")
    plt.xlabel("Time (ns)")
    plt.title(
        "Actual and ideal values of slopes of cumulative absolute velocity distributions"
    )
    plt.legend()
    plt.show()
    # Calculate mean and standard deviations for plotting a normal distribution overlay
    mean_position = np.mean(particle_positions)
    position_std_dev = np.std(particle_positions)
    # print("Mean", mean_position)
    # print("Standard deviation", position_std_dev)
    # Create a list of all position values
    position_values = np.linspace(0, 1000)
    # Call the normal distribution on all position values
    position_normal_distribution = [
        normal_distribution(any_value, mean_position, position_std_dev)
        for any_value in position_values
    ]
    # print(position_values)
    # print(position_normal_distribution)
    # Plotting Histogram of particle_positions
    plt.hist(particle_positions, bins=50, label="Positions")
    plt.plot(
        position_values,
        position_normal_distribution,
        label="Normal distribution overlay",
    )
    plt.ylabel("Frequency")
    plt.xlabel("Position(nm)")
    plt.title(
        f"Histogram for particle position in Brownian Motion on a {line_length} nm Line"
    )
    plt.legend()
    plt.show()

    # Plotting Histogram of particle_velocities
    plt.hist(particle_velocity, label="Velocity", color="orange")
    plt.ylabel("Frequency")
    plt.xlabel("Velocity (nm/s)")
    plt.title(
        f"Histogram for particle velocity in Brownian Motion on a {line_length} nm Line"
    )
    plt.legend()
    plt.show()
    position_slope = []
    slope_time = []  # to store all the times the slope is calculated for
    # loop to calculate slope of the position with time
    for i in range(0, len(particle_positions) - 1, 100):
        position_change = particle_positions[i + 1] - particle_positions[i]
        time_change = particle_time[i + 1] - particle_time[i]
        slope = position_change / time_change
        position_slope.append(slope)
        slope_time.append(particle_time[i])
    # Plotting the slope of position graph
    plt.scatter(slope_time, position_slope, label="Position Slope")
    plt.ylabel("Position Slope (nm/ns)")
    plt.xlabel("Time (ns)")
    plt.legend()
    plt.title(
        f"Slope of Position vs Time for Brownian Motion on a {line_length} nm Line"
    )

    plt.show()


def gaussian(last_position, sigma):
    # generate a value in a gaussian distribution based on the mean being the previous position/velocity
    x = np.random.normal(last_position, sigma)
    return x


def normal_distribution(current_position, mean, std_dev):
    return (
        5
        * 10**4
        * 1
        / (math.sqrt(2 * math.pi) * std_dev)
        * math.exp(-0.5 * math.pow(((current_position - mean) / std_dev), 2))
    )


def calculate_slope(x1, y1, x2, y2):
    return (y2 - y1) / (x2 - x1)


if __name__ == "__main__":
    brownian_motion_simulation()
