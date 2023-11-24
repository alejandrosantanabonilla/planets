import numpy as np
import os
from tqdm import tqdm

def objective_function(x, y):
    """Objective function to be minimized"""
    return (x - 3.14)**2 + (y - 2.72)**2 + np.sin(3 * x + 1.41) + np.sin(4 * y - 1.73)
   
    
def plot_contour(x_min, x_max, y_min, y_max, coordinates, idx):
    """Plot the objective function as a contour plot.

    Args:
        x_min: Minimum value of the x-axis
        x_max: Maximum value of the x-axis
        y_min: Minimum value of the y-axis
        y_max: Maximum value of the y-axis
        coordinates: List of coordinates of the current best solution
        idx: Index of the current best solution

    Returns:
        fig, ax: The figure and axis of the contour plot
    """
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    x, y = np.meshgrid(np.linspace(x_min, x_max, 1000), np.linspace(y_min, y_max, 1000))
    z = objective_function(x, y)

    fig, ax = plt.subplots()
    ax.contourf(x, y, z, cmap='viridis')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Objective Function')
    fig.colorbar(ax.collections[0], label='Objective Value')

    ax.scatter(coordinates[0], coordinates[1], color='red', marker='o', label=''.join(["Coordinate_",str(idx)]))
    ax.legend()
    
    return fig, ax

def update_particles(n_particles, c1, c2, w, number_gen):
    """Update particle positions and velocities

    Args:
        n_particles: Number of particles for the swarm
        c1: First coeffient of the PSO velocity
        c2: Second coefficient for the PSO velocity
        w:  Inertia variable
        number_gen: Number of generations

    Returns:
        gbest: List of all the gbest values throughout the generations

    """
    import json
    
    gbest_hist=[]

    V = np.random.randn(2, n_particles) * 0.1
    X = np.random.rand(2, n_particles) * 10

    pbest = X.copy()
    pbest_obj = np.zeros(n_particles)
    
    for i in range(n_particles):
        pbest_obj[i] = objective_function(X[0, i], X[1, i])

    gbest = pbest[:, pbest_obj.argmin()]
    gbest_obj = pbest_obj.min()

    # Optimization loop
    for iteration in tqdm(range(int(number_gen)),desc="PSO iter: "):

        # Update particle velocities
        r1, r2 = np.random.rand(2)
        V = w * V + c1 * r1 * (pbest - X) + c2 * r2 * (gbest.reshape(-1, 1) - X)

        # Update particle positions
        X = X + V

        # Calculate objective function values for current positions
        obj = objective_function(X[0], X[1])

        # Update personal best positions and objective values
        pbest[:, (pbest_obj >= obj)] = X[:, (pbest_obj >= obj)]
        pbest_obj = np.array([pbest_obj, obj]).min(axis=0)

        # Update global best position and objective value
        gbest = pbest[:, pbest_obj.argmin()]
        gbest_obj = pbest_obj.min()

        # Save X and V variables to a JSON file for subsequent plotting
        with open(''.join(["particle_",str(iteration),".json"]), 'w') as f:
            json.dump({'X': X.tolist(), 'V': V.tolist()}, f)

        #Print positions and velocities
        #print("Iteration:", iteration + 1)
        #print("Particle Positions:")

        #for i in range(n_particles):
        #    print("Particle", i + 1, ":", X[:, i])

        #print("\nParticle Velocities:")
        #for i in range(n_particles):
        #    print("Particle", i + 1, ":", V[:, i])

        #print("\nGlobal Best Position:", gbest)
        gbest_hist.append(gbest)
      
        print("Global Best Objective Value:", gbest_obj)
        print("\n")

    return gbest_hist

def s_m_files(filename, extension):
    """This function moves all files with the specified extension from the current directory 
       to a temporary directory with the same name as the filename.

    Args:
       filename: The name of the file that will be used to create the temporary directory
       extension: The file extension of the files to be moved

    Returns:
       None
    """
    import os
    import shutil

    # Get the current directory
    current_dir = os.getcwd()

    # Create a temporary directory
    temp_dir = os.path.join(current_dir, str(filename))
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # Find all .json files in the current directory
    for filename in tqdm(os.listdir(current_dir),desc="Moving files: "):
        if filename.endswith(''.join([".",str(extension)])):
          # Move the .json file to the temporary directory
            shutil.move(os.path.join(current_dir, filename), os.path.join(temp_dir, filename))

    print('All {} files have been moved to the temporary directory: {}'.format(extension, temp_dir))

def create_gif(png_files, gif_filename):
    """Creates a GIF animation from a list of PNG files.

    Args:
        png_files: List of PNG file paths.
        gif_filename: Path to save the GIF animation.
    """
    from PIL import Image, ImageDraw

    img, *imgs  = [Image.open(fn) for fn in png_files]
    img.save(gif_filename, format='GIF', append_images=imgs,
             save_all=True, loop=0)
    

def pso_simulation(n_particles, c1, c2, w, number_gen, x_min, y_min, x_max, y_max, image_name="contour_plot", animation_name="PSO"):
    """ Particle Swarm Optimisation evolution """

    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    coordinates = update_particles(n_particles, c1, c2, w, number_gen)

    for idx, coord in enumerate(coordinates):
      fig, ax = plot_contour(float(x_min), float(x_max), float(y_min), float(y_max), coord, idx)
      fig.savefig(''.join([str(image_name),str(idx),'.png']))
      plt.close(fig)

    gif_filename = ''.join([str(animation_name),'.gif'])
    png_files = [''.join([str(image_name),str(idx),'.png']) for idx in range(0,len(coordinates))]
    create_gif(png_files, gif_filename)

    s_m_files("temp_json","json")
    s_m_files("temp_png", "png")
    
if __name__ == "__main__":

    # Simulation parameters
    n_particles = 15
    c1 = c2 = 0.1
    w = 0.8
    number_gen=15
    
    x_min, x_max = 0, 10
    y_min, y_max = 0, 10

    pso_simulation(n_particles, c1, c2, w, number_gen, x_min, y_min, x_max, y_max, image_name="contour_plot", animation_name="PSO")
