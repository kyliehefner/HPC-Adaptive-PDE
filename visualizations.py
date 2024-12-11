import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# folder containing the output files
data_folder = "HPC-Adaptive-PDE/data"

# list of time steps to visualize
time_steps = range(0, 200, 10)
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
    plt.imshow(grid, cmap="hot", interpolation="nearest", vmin=0.0, vmax=1.0)
    plt.colorbar(label="Concentration", shrink=0.7)
    plt.title(f"Step {step}")

plt.tight_layout()
plt.show()


# Plot contour plots for each time step
plt.figure(figsize=(15, num_rows * 3))
for idx, step in enumerate(time_steps):
    file_path = os.path.join(data_folder, f"output_step_{step}.csv")
    grid = np.loadtxt(file_path, delimiter=",")
    
    plt.subplot(num_rows, num_cols, idx + 1)  # Adjust rows/columns based on the number of time steps
    plt.contourf(grid, cmap="viridis", levels=20)
    plt.colorbar(label="Concentration", shrink=0.8)
    plt.title(f"Step {step}")

plt.tight_layout()
plt.savefig("contour_evolution.png", dpi=300)
plt.show()


# Initialize the figure
fig, ax = plt.subplots(figsize=(8, 8))
cax = ax.imshow(np.zeros((100, 100)), cmap="hot", vmin=0.0, vmax=1.0)
plt.colorbar(cax, label="Concentration")
ax.set_title("Advection-Diffusion Evolution")

# Update function for animation
def update(frame):
    grid = np.loadtxt(os.path.join(data_folder, f"output_step_{frame}.csv"), delimiter=",")
    cax.set_data(grid)
    ax.set_title(f"Step {frame}")
    return cax,

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=time_steps, interval=200)
ani.save("solution_evolution.gif", writer='pillow', fps=5)  # Save as GIF
plt.show()

