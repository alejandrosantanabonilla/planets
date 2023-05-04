import glob
import os
from operations import *
import shutil

def read_input_file(input_file):
    """ Function to read a valid input file for searching 
        structures.
    """

    with open(str(input_file)) as f:
       for line in f:
          if line.startswith("startstep"):
             temp_st = line.partition('=')
         
          if line.startswith("maxstep"):
             temp_mt = line.partition('=')
         
          if line.startswith("popsize"):
             temp_ps = line.partition('=')

          if line.startswith("command"):
             temp_com = line.partition('=')

    return int(temp_st[-1]), int(temp_mt[-1]), int(temp_ps[-1]), str(temp_com[-1])

def make_config_folders(gen_dir, parent_dir, number, name_new_str):
    """ Function to create folders containing different 
        configurations
    """
    
    conf_dir="config_{}".format(number)

    path = os.path.join(gen_dir, conf_dir)

    try: 
        os.mkdir(path) 

    except OSError as error: 
        print(error)

    dst_path=os.path.join(conf_dir,path)
    src_path=os.path.join(gen_dir, name_new_str)
    shutil.move(src_path, dst_path)

def sub_list(gen_dir):
  """ Function to
  """
  sub=[f.path for f in os.scandir(gen_dir) if f.is_dir()]

  return sub

def vasp_energy(command, energy, subfolder_list, output_name="pop_res.json"):
    """ Function to perform VASP calculations
        and obtain the energies. It is a prototype
    """
    import json
    
    keys=[]; eng=[]
    
    for idx, values in enumerate(subfolder_list):
        s = np.random.uniform(0,1)
        
        completeName=os.path.join(values,"OUTCAR")

        #Fake energy from random number
        with open(completeName,"w") as file1:
            file1.write("{}".format(energy+s))

        #Name of the keys based on generation 
        keys.append('_'.join(map(lambda x: x, values.split(os.sep)[-2:])))

        with open(completeName, "r") as f:
            eng.append([f.readlines()])


    #Path to store the dictionary gen_1
    head_tail=os.path.split(subfolder_list[0])
    fullname=os.path.join(head_tail[0],output_name)

    #Cleaning and creating the dictionary
    last=sum(sum(eng,[]),[])
    dictionary=dict(zip(keys,last))

    with open(fullname, "w") as outfile:
       json.dump(dictionary, outfile)

      
def main():

    startstep, maxstep, pop, command=read_input_file("input.dat")

    #Read the cluster.POSCAR file and convert to Cartesian with selective dynamics style
    clusters = read("cluster.poscar", format='vasp')
    substrate = read("substrate.poscar", format='vasp')
   
    parent_dir = os.getcwd()
    gen_count=startstep

    values=global_random_numbers(pop, maxstep, seed=15)
    velocities=random_velocities(pop, maxstep, seed=42)
    z_min=2.5    
    
    while gen_count < maxstep:
        
        directory="gen_{}".format(str(gen_count))
        path = os.path.join(parent_dir, directory)
        
        try: 
            os.mkdir(path)
            
        except OSError as error: 
            print(error)
        
        if gen_count <= 1:
            
           os.chdir(path)
           gen_dir=os.getcwd()
           
           for i in range(pop):
             name_new_str="".join(["combined_",str(i),".poscar"])
             structure_creator(clusters, substrate, name_new_str, values[gen_count-1][i], z_min)
             make_config_folders(gen_dir, parent_dir, i, name_new_str)

           subfolders=sub_list(gen_dir)
           vasp_energy("blah", -187.9896, subfolders)
           os.chdir(parent_dir)
             
        else:
           os.chdir(path)
           gen_dir_new=os.getcwd()

           for i in range(pop):
             name_new_str_rest="".join(["combined_",str(i),".poscar"])
             structure_creator(clusters, substrate, name_new_str_rest, values[gen_count][i], z_min)
             make_config_folders(gen_dir_new, parent_dir, i, name_new_str_rest)

           subfolders_new=sub_list(gen_dir_new)
           
           os.chdir(parent_dir)

        gen_count += 1


if __name__ == "__main__":
    
   main()
