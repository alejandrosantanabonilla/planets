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

├── cluster.poscar
├── gen_1
│   ├── config_0
│   │   ├── combined_0.poscar
│   │   └── OUTCAR
│   ├── config_1
│   │   ├── combined_1.poscar
│   │   └── OUTCAR
│   ├── config_2
│   │   ├── combined_2.poscar
│   │   └── OUTCAR
│   ├── config_3
│   │   ├── combined_3.poscar
│   │   └── OUTCAR
│   ├── config_4
│   │   ├── combined_4.poscar
│   │   └── OUTCAR
│   ├── config_5
│   │   ├── combined_5.poscar
│   │   └── OUTCAR
│   ├── config_6
│   │   ├── combined_6.poscar
│   │   └── OUTCAR
│   ├── config_7
│   │   ├── combined_7.poscar
│   │   └── OUTCAR
│   ├── config_8
│   │   ├── combined_8.poscar
│   │   └── OUTCAR
│   ├── config_9
│   │   ├── combined_9.poscar
│   │   └── OUTCAR
│   └── pop_res.json
├── gen_2
│   ├── config_0
│   │   └── combined_0.poscar
│   ├── config_1
│   │   └── combined_1.poscar
│   ├── config_2
│   │   └── combined_2.poscar
│   ├── config_3
│   │   └── combined_3.poscar
│   ├── config_4
│   │   └── combined_4.poscar
│   ├── config_5
│   │   └── combined_5.poscar
│   ├── config_6
│   │   └── combined_6.poscar
│   ├── config_7
│   │   └── combined_7.poscar
│   ├── config_8
│   │   └── combined_8.poscar
│   └── config_9
│       └── combined_9.poscar
├── gen_3
│   ├── config_0
│   │   └── combined_0.poscar
│   ├── config_1
│   │   └── combined_1.poscar
│   ├── config_2
│   │   └── combined_2.poscar
│   ├── config_3
│   │   └── combined_3.poscar
│   ├── config_4
│   │   └── combined_4.poscar
│   ├── config_5
│   │   └── combined_5.poscar
│   ├── config_6
│   │   └── combined_6.poscar
│   ├── config_7
│   │   └── combined_7.poscar
│   ├── config_8
│   │   └── combined_8.poscar
│   └── config_9
│       └── combined_9.poscar
├── gen_4
│   ├── config_0
│   │   └── combined_0.poscar
│   ├── config_1
│   │   └── combined_1.poscar
│   ├── config_2
│   │   └── combined_2.poscar
│   ├── config_3
│   │   └── combined_3.poscar
│   ├── config_4
│   │   └── combined_4.poscar
│   ├── config_5
│   │   └── combined_5.poscar
│   ├── config_6
│   │   └── combined_6.poscar
│   ├── config_7
│   │   └── combined_7.poscar
│   ├── config_8
│   │   └── combined_8.poscar
│   └── config_9
│       └── combined_9.poscar
├── gen_5
│   ├── config_0
│   │   └── combined_0.poscar
│   ├── config_1
│   │   └── combined_1.poscar
│   ├── config_2
│   │   └── combined_2.poscar
│   ├── config_3
│   │   └── combined_3.poscar
│   ├── config_4
│   │   └── combined_4.poscar
│   ├── config_5
│   │   └── combined_5.poscar
│   ├── config_6
│   │   └── combined_6.poscar
│   ├── config_7
│   │   └── combined_7.poscar
│   ├── config_8
│   │   └── combined_8.poscar
│   └── config_9
│       └── combined_9.poscar
├── gen_6
│   ├── config_0
│   │   └── combined_0.poscar
│   ├── config_1
│   │   └── combined_1.poscar
│   ├── config_2
│   │   └── combined_2.poscar
│   ├── config_3
│   │   └── combined_3.poscar
│   ├── config_4
│   │   └── combined_4.poscar
│   ├── config_5
│   │   └── combined_5.poscar
│   ├── config_6
│   │   └── combined_6.poscar
│   ├── config_7
│   │   └── combined_7.poscar
│   ├── config_8
│   │   └── combined_8.poscar
│   └── config_9
│       └── combined_9.poscar
├── gen_7
│   ├── config_0
│   │   └── combined_0.poscar
│   ├── config_1
│   │   └── combined_1.poscar
│   ├── config_2
│   │   └── combined_2.poscar
│   ├── config_3
│   │   └── combined_3.poscar
│   ├── config_4
│   │   └── combined_4.poscar
│   ├── config_5
│   │   └── combined_5.poscar
│   ├── config_6
│   │   └── combined_6.poscar
│   ├── config_7
│   │   └── combined_7.poscar
│   ├── config_8
│   │   └── combined_8.poscar
│   └── config_9
│       └── combined_9.poscar
├── gen_8
│   ├── config_0
│   │   └── combined_0.poscar
│   ├── config_1
│   │   └── combined_1.poscar
│   ├── config_2
│   │   └── combined_2.poscar
│   ├── config_3
│   │   └── combined_3.poscar
│   ├── config_4
│   │   └── combined_4.poscar
│   ├── config_5
│   │   └── combined_5.poscar
│   ├── config_6
│   │   └── combined_6.poscar
│   ├── config_7
│   │   └── combined_7.poscar
│   ├── config_8
│   │   └── combined_8.poscar
│   └── config_9
│       └── combined_9.poscar
├── gen_9
│   ├── config_0
│   │   └── combined_0.poscar
│   ├── config_1
│   │   └── combined_1.poscar
│   ├── config_2
│   │   └── combined_2.poscar
│   ├── config_3
│   │   └── combined_3.poscar
│   ├── config_4
│   │   └── combined_4.poscar
│   ├── config_5
│   │   └── combined_5.poscar
│   ├── config_6
│   │   └── combined_6.poscar
│   ├── config_7
│   │   └── combined_7.poscar
│   ├── config_8
│   │   └── combined_8.poscar
│   └── config_9
│       └── combined_9.poscar
├── input.dat
├── operations.py





