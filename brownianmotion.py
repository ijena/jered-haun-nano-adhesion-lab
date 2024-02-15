import matplotlib.pyplot as plt
import math


# TO DO add error handling for inputs
# def value_error_handling(user_input, data_type):
#    try:

line_length = (float)(
    input("Enter length of the box in nanometre (nm)\n")
)  # input from user determines length of the box

x_position = line_length / 2  # initial position of the particle at midpoint of the line
total_time = (int)(
    input("Enter total time of simulation in nanoseconds (ns)\n")
)  # input from user determines time of simulation
particle_radius = (float)(input("Enter radius of the particle in nanometre\n"))
particle_mass = (float)(input("Enter mass of particle\n"))
# lower_bound and upper_bound refer to the endpoints of the line
lower_bound = 0
upper_bound = line_length
particle_positions = [
    x_position
]  # particle_positions is created to store all the positions of the particle so that it can be plot later. line_length/2 nm is the initial position of the particle at time = 0 ns
particle_time = [
    0
]  # particle_time is created to store the time passed since start of simulation. 0ns is the start of the simulation

# particle radius is 100 nm
viscosity_liquid = 0.01  # 0.01 poise is viscosity of water
inverse_viscous_relaxation_time = (
    3 * math.pi * 2 * particle_radius
) / particle_mass  # calculating inverse viscous relaxation time using equation 4 from Hammer- English paper
for time in range(1, total_time):
    # use equation from research paper to determine motion
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
    particle_time.append(time)

# plot particle position vs time
plt.plot(particle_time, particle_positions)
plt.ylabel("Position of particle in nanometre (nm)")
plt.xlabel("Time passed since start of simulation in nanoseconds (ns)")
plt.title(f"Simulation of a Brownian motion on a {line_length} nanometre line")
plt.show()
