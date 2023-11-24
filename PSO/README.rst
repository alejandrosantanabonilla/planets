**Tutorial to use the Particle Swarm Optimization python implemention:**

**Introduction**

Particle Swarm Optimization (PSO) is a population-based metaheuristic optimization algorithm inspired by the social behavior of birds flocking, fish schooling, or insect swarming. 
PSO is a robust and efficient optimization technique that has been successfully applied to a wide range of problems in various fields, including engineering, science, and business.

This tutorial will provide a step-by-step guide on how to use the PSO implementation in Python.

**Prerequisites**

* Basic understanding of Python programming
* Familiarity with optimization algorithms

**Installation Process:**

Before proceeding with the tutorial, ensure that you have installed the necessary libraries, including:

* NumPy
* Matplotlib
* tqdm
* json
* PIL

**Installing Dependencies:**

To insall libraries such as NumPy, JSON, tqdm, and PIL please use the following commands:

.. code-block:: bash

   pip install numpy
   pip install json   
   pip install matplotlib
   pip install json
   pip install Pillow
   pip install tqdm

**Note:** Make sure you have pip installed on your system. If not, you can install it using the following command:

.. code-block:: bash

   sudo apt update
   sudo apt install python3-pip

**How to run:**

An example can be run directly after the code has been downloaded and libraries have been installed. This is done using
this command:

.. code-block:: bash

  python pso_code.py

after running this command, 2 folders (temp_json, temp_png) are going to be created. You can open the *.json* file with any
text editor (emacs, nano, vim, etc) and see the positions and velocities of the particles per each iteration. The **PSO.gif**
file contains a small animation showing the best position (gbest) found by the algorithm at that iteration. Values will be
printed in the screen. The animation can be seen by using the following command:

.. code-block:: bash

   firefox PSO.gif

**What to change:**

To change the parameters of the code, you can find them in the section:

.. code-block:: python

   if __name__ == "__main__":

       # Simulation parameters
       n_particles = 25
       c1 = c2 = 0.1
       w = 0.8
       number_gen=15

       #User defined grid for n_particle search
       x_min, x_max = 0, 10
       y_min, y_max = 0, 10

       #Running the PSO simulation
       pso_simulation(n_particles, c1, c2, w, number_gen, x_min, y_min, x_max, y_max, image_name="contour_plot", animation_name="PSO")

The PSO name convention is used throughout the code (omega, c1, c2, number of particles, and number of generations). The rest of the variables
are the dimensions of the grid for the objective function ((xmin,xmax) X (ymin,ymax)). To change the objective function, you just need to change
the pytho function:

.. code-block:: python

    def objective_function(x, y):
    """Objective function to be minimized"""
    return (x - 3.14)**2 + (y - 2.72)**2 + np.sin(3 * x + 1.41) + np.sin(4 * y - 1.73)

**Conclusion**

This tutorial provides a basic implementation of PSO in Python. For more advanced applications, you can explore 
various extensions and modifications of the algorithm in terms of new inertia functions, velocity clamping 
among others.

**Next steps:**

# Step 1. blala
                 
# Step 2. blala


