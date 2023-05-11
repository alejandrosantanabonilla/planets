import os
import subprocess as sp
from subprocess import call

num_cores = int(1) if num_cores is None else int(num_cores)
init_dir=os.getcwd()
    
out_file = os.path.join(init_dir,"output.log")

cmd = (
        f'mpirun -np {num_cores}'
        f'{vasp_path}'
       )

#Note that sp.call will hold the program until completion
# of the calculation.

try:
   with open(out_file, 'w') as result:
       sp.call(
             cmd,
             stdin=sp.PIPE,
             stdout=result,
             stderr=sp.PIPE,
             shell=True,
             cwd=init_dir
           )
finally:
        pass
