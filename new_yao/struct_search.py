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

if __name__ == "__main__":
    startstep,maxstep,pop,command=read_input_file("input.dat")

    gen_count=startstep
    while gen_count < maxstep :
        print ("gen_{}".format(str(gen_count)))
        gen_count += 1
