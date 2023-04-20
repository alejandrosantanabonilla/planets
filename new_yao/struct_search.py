import glob
import os


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

def make_config_folders(pop, gen_dir, parent_dir):
    """ Function to create folders containing different 
        configurations
    """
    
    for i in range(pop):
       conf_dir="config_{}".format(i)
       path = os.path.join(gen_dir, conf_dir)

       try: 
           os.mkdir(path) 
       except OSError as error: 
           print(error)

    os.chdir(parent_dir)
    print ("Succesfully created {} configuration folders".format(pop))


def random_vel_list(step=0.3):
    """ Function to create a random list
    """

    import random
    
    random_v1 = random.uniform(-step, step)
    random_v2 = random.uniform(-step, step)
    random_v3 = random.uniform(-step, step)
    random_v4 = random.uniform(-step, step)
    random_v5 = random.uniform(-step, step)

    with open('velocity.dat', 'w') as f:
        f.write(f'{random_v1}\n')
        f.write(f'{random_v2}\n')
        f.write(f'{random_v3}\n')
        f.write(f'{random_v4}\n')
        f.write(f'{random_v5}\n')

    print('velocity\n',random_v1, random_v2, random_v3, random_v4, random_v5)

def random_rot_list(step=1.0):
    """ Function to obtain random numbers to perform 
        an arbitrary rotation.
    """
    import random
    
    random_a = random.uniform(-step, step)
    random_b = random.uniform(-(step/2.0), step/2.0)
    random_c = random.uniform(-step, step)
    random_x = random.random()
    random_y = random.random()

    with open('structure.param', 'w') as f:
        f.write(f'{random_a}\n')
        f.write(f'{random_b}\n')
        f.write(f'{random_c}\n')
        f.write(f'{random_x}\n')
        f.write(f'{random_y}\n')

    print('random for rotation\n',random_a, random_b, random_c,'\nrandom for displacement\n', random_x, random_y)
    
if __name__ == "__main__":
    
    startstep, maxstep, pop, command=read_input_file("input.dat")

    parent_dir = os.getcwd()
    gen_count=startstep
    
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
           make_config_folders(pop, gen_dir, parent_dir)

           # Calculation
           random_vel_list()
           random_rot_list()
           
           
        else:
           os.chdir(path)
           gen_dir_new=os.getcwd()

           make_config_folders(pop, gen_dir_new, parent_dir)

           
        gen_count += 1
