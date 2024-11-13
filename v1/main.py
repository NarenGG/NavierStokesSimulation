import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import cv2
import glob
import re

# Initialize the pressure grid with higher values on the left side
definition = 50
pressure_grid = np.zeros((definition + 2, definition + 2, definition + 2))
pressure_grid[:, :, :definition//2 + 1] = 1.0
velocity_grid = np.zeros((3, definition + 2, definition + 2, definition + 2))
rho = 1.0
mu = 0.01
dt = 0.1
dx = dy = dz = 1.0
num_time_steps = 1000
pressure_grids, velocity_grids = [pressure_grid], [velocity_grid]

def calculate_next_state(pressure_grid, velocity_grid, rho, mu, dt, dx, dy, dz):
    """
    Calculate the next state of the pressure and velocity grids using the Navier-Stokes equation.
    
    Parameters:
    pressure_grid (numpy.ndarray): 3D array of pressures.
    velocity_grid (numpy.ndarray): 4D array of velocities (shape: (3, nx, ny, nz)).
    rho (float): Density of the fluid.
    mu (float): Dynamic viscosity of the fluid.
    dt (float): Time step.
    dx (float): Grid spacing in x direction.
    dy (float): Grid spacing in y direction.
    dz (float): Grid spacing in z direction.
    
    Returns:
    tuple: Updated pressure grid and velocity grid.
    """
    nx, ny, nz = pressure_grid.shape
    u, v, w = velocity_grid
    
    # Initialize the next state pressure and velocity grids
    next_pressure_grid = np.copy(pressure_grid)
    next_velocity_grid = np.copy(velocity_grid)
    
    # Replace NaN values with 0
    pressure_grid = np.nan_to_num(pressure_grid, nan=0.0)
    u = np.nan_to_num(u, nan=0.0)
    v = np.nan_to_num(v, nan=0.0)
    w = np.nan_to_num(w, nan=0.0)
    
    # Define the maximum and minimum values to clip to
    max_value = np.finfo(np.float64).max
    min_value = np.finfo(np.float64).min
    
    # Calculate the gradients and Laplacian
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                # Calculate the velocity gradients
                du_dx = (u[i+1, j, k] - u[i-1, j, k]) / (2 * dx)
                dv_dy = (v[i, j+1, k] - v[i, j-1, k]) / (2 * dy)
                dw_dz = (w[i, j, k+1] - w[i, j, k-1]) / (2 * dz)
                
                # Calculate the Laplacian of pressure
                laplacian_p = (
                    (pressure_grid[i+1, j, k] - 2 * pressure_grid[i, j, k] + pressure_grid[i-1, j, k]) / dx**2 +
                    (pressure_grid[i, j+1, k] - 2 * pressure_grid[i, j, k] + pressure_grid[i, j-1, k]) / dy**2 +
                    (pressure_grid[i, j, k+1] - 2 * pressure_grid[i, j, k] + pressure_grid[i, j, k-1]) / dz**2
                )
                
                # Update the pressure using the Navier-Stokes equation
                next_pressure_grid[i, j, k] = pressure_grid[i, j, k] + dt * (
                    -rho * (u[i, j, k] * du_dx + v[i, j, k] * dv_dy + w[i, j, k] * dw_dz) +
                    mu * laplacian_p
                )
                
                # Update the velocity components
                next_velocity_grid[0, i, j, k] = np.clip(u[i, j, k] + dt * (
                    -u[i, j, k] * du_dx - v[i, j, k] * dv_dy - w[i, j, k] * dw_dz +
                    mu * (laplacian_p / rho)
                ), min_value, max_value)
                
                next_velocity_grid[1, i, j, k] = np.clip(v[i, j, k] + dt * (
                    -u[i, j, k] * du_dx - v[i, j, k] * dv_dy - w[i, j, k] * dw_dz +
                    mu * (laplacian_p / rho)
                ), min_value, max_value)
                
                next_velocity_grid[2, i, j, k] = np.clip(w[i, j, k] + dt * (
                    -u[i, j, k] * du_dx - v[i, j, k] * dv_dy - w[i, j, k] * dw_dz +
                    mu * (laplacian_p / rho)
                ), min_value, max_value)
    
    return next_pressure_grid, next_velocity_grid

def plot_3d_grid(grid, title):
    """
    Plot a 3D grid using Matplotlib.
    
    Parameters:
    grid (numpy.ndarray): 3D array to plot.
    title (str): Title of the plot.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Create a grid of coordinates
    x, y, z = np.indices(grid.shape)
    
    # Flatten the arrays
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()
    values = grid.flatten()
    
    # Create a scatter plot
    scatter = ax.scatter(x, y, z, c=values, cmap='viridis')
    
    # Add color bar
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=5)
    cbar.set_label('Pressure')
    
    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    plt.show()

def plot_slices(grid, title):
    """
    Plot slices of a 3D grid using Matplotlib.
    
    Parameters:
    grid (numpy.ndarray): 3D array to plot.
    title (str): Title of the plot.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    slices = [grid[grid.shape[0] // 2, 1:-1, 1:-1], grid[1:-1, grid.shape[1] // 2, 1:-1], grid[1:-1, 1:-1, grid.shape[2] // 2]]
    slice_titles = ['XY Slice', 'XZ Slice', 'YZ Slice']
    
    for ax, slice, slice_title in zip(axes, slices, slice_titles):
        cax = ax.imshow(slice, cmap='viridis', origin='lower', vmin=0.0, vmax=1.0)
        ax.set_title(slice_title)
        fig.colorbar(cax, ax=ax)
    
    fig.suptitle(title)
    plt.savefig(title)
    plt.close()

# Create a directory to store the data
os.system('rm -rf ./data && mkdir data')

# Plot the initial state
plot_slices(pressure_grids[-1], f"data/0")

# Run the simulation for the specified number of time steps
for epoch in range(1, num_time_steps + 1):
    next_pressure_grid, next_velocity_grid = calculate_next_state(pressure_grids[-1], velocity_grids[-1], rho, mu, dt, dx, dy, dz)
    next_pressure_grid = np.nan_to_num(next_pressure_grid, nan=0)
    next_velocity_grid = np.nan_to_num(next_velocity_grid, nan=0)
    pressure_grids.append(next_pressure_grid)
    velocity_grids.append(next_velocity_grid)
    print(f"Processed Grid at Epoch {epoch}")
    plot_slices(pressure_grids[-1], f"data/{epoch}")

def extract_number(filename):
    """
    Extract the number from the filename.
    
    Parameters:
    filename (str): Filename to extract the number from.
    
    Returns:
    int: Extracted number.
    """
    match = re.search(r'(\d+)\.png$', filename)
    return int(match.group(1)) if match else float('inf')

# Get the list of filenames and sort them by the number before .png
filenames = sorted(glob.glob('data/*.png'), key=extract_number)

# Read the images and store them in an array
img_array = []
for filename in filenames:
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width, height)
    img_array.append(img)

# Create a video from the images
out = cv2.VideoWriter('output_video.avi', cv2.VideoWriter_fourcc(*'DIVX'), 15, size)
for img in img_array:
    out.write(img)
out.release()
