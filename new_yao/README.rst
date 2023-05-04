Instructions
==================


To use ase, you just need to install using pip in your computer. This is done by

.. code-block:: python

   pip install ase
   
Please download all the files here to use both codes. 


How to use the code?
=========================

To use the code the files **struc_search.py** and **operations.py** should be place in the same folder. To get the creation of the 
folders and the different structures organised in different folders, you need to type the following:

.. code-block:: python

   python struc_search.py
   
This will produce the following set of folders:

.. code-block:: bash

   ── cluster.poscar
   ├── gen_1
   │   ├── config_0
   │   │   ├── combined_0.poscar
   │   │   └── OUTCAR
   │   ├── config_1
   │   │   ├── combined_1.poscar
   │   │   └── OUTCAR
            .
            .
            .
            .
   └── pop_res.json
   ├── gen_2
            .
            .
            .
    ├── gen_9
   

This is for **9 generations** and in each generation we have **10 individuals for a population**.

There are two special files which are called **vel_param.npy** and **struc_param.npy**. They contain the randomly uniform distributed parameters
from numpy. Unfortunately, I have found out that calling this function every single cycle produces a correlation and therefore you are not 
creating "truly" random numbers. We need to compute them before the cycle starts and then take from those two matrices. They are in **binary**
format.

Workflow
===============

The workflow now is completed. I can create the structures from random numbers allocate them into folders and create a **fake energy" whose 
values and the name of the folder of origin are stored in a json file named "pop_res.json". The info can be seen using 


.. code-block:: bash
 
   head pop_res.json







