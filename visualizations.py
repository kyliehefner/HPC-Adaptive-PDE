import numpy as np
import matplotlib.pyplot as plt
import os

# folder containing the output files
data_folder = "HPC-Adaptive-PDE/data"

# list of time steps to visualize
time_steps = range(0, 100, 10)
num_steps = len(time_steps)

# determine the grid layout for subplots (rows and columns)
num_cols = 5
num_rows = (num_steps + num_cols - 1) // num_cols  # compute rows dynamically

# create the figure for subplots
plt.figure(figsize=(15, num_rows * 3))

# plot the solution at each time step
for idx, step in enumerate(time_steps):
    file_path = os.path.join(data_folder, f"output_step_{step}.csv")
    grid = np.loadtxt(file_path, delimiter=",")
    
    # add a subplot for each time step
    plt.subplot(num_rows, num_cols, idx + 1)  # use dynamic grid positioning
    plt.imshow(grid, cmap="hot", interpolation="nearest")
    plt.colorbar(label="Concentration", shrink=0.7)
    plt.title(f"Step {step}")

plt.tight_layout()
plt.show()

