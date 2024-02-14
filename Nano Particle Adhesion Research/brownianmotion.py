import random
import matplotlib.pyplot as plt

# change all values that get affect by line_length
line_length = 1000  # length of the line in nanometre (nm)
x_position = line_length / 2  # initial position of the particle at midpoint of the line
total_time = 1000  # this simulation will run for 1000 seconds
# lower_bound and upper_bound refer to the endpoints of the line
lower_bound = 0
upper_bound = 500
# particle_positions is created to store all the positions of the particle so that it can be plot later. 250 nm is the position of the particle at t = 0s
particle_positions = [250]
# particle_time is created to store the time passed since start of simulation. 0s is the start of the simulation
particle_time = [0]
# step_size refers to how far can the particle move in 1 second
step_size = 1
# calculate in nanoseconds instead of seconds
# there is no fixed position step_size

# particle radius is 100 nm
for time in range(1, total_time):
    # randomly chooses whether the particle moves backwards or forwards in that second
    x_position = random.choice(
        [x_position - step_size, x_position + step_size]
    )  # change to the formula from the research paper
    # checks if the particle remains within the bounds
    if x_position < lower_bound:
        x_position += 2 * step_size
    elif x_position > upper_bound:
        x_position -= 2 * step_size

    particle_positions.append(x_position)
    particle_time.append(time)

# plot particle position vs time
plt.plot(particle_time, particle_positions)
plt.ylabel("Position of particle in nanometre (nm)")
plt.xlabel("Time passed since start of simulation in seconds (s)")
plt.title("Simulation of a Brownian motion on a 500 nanometre line")
plt.show()
