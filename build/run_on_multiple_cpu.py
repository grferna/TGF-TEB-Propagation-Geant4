#!/usr/bin/env python3

import platform
from subprocess import call
from mpi4py import MPI
from functools import partial
from multiprocessing.dummy import Pool
import random
import datetime

## Initialization of MPI routines
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
##

print('Number of threads '+ str(nprocs))

# number of CPU to use (for local run)
NB_CPU_LOCAL = 5

# number of times the list of command is copied (usually usefull for large number of CPU, like for clusters)
NB_copy_commands = 100

nr_record_to_get_per_run = 1000                 

# lists of parameters to test
initial_alt_list = [18.,16.,14.,12.,10.]
beaming_angle_list = [5., 10.,20.,30.,40.,50.,60.]
# tilt_angle_list = [5.,10.,20.,30.,40.,50.,60.]
tilt_angle_list = [0.]
#beaming_type_list=["Gaussian","Uniform"]
beaming_type_list=["Uniform"]
source_sigma_time_list=[0.0]

# defining a list of commands to be run in parallel

commands=[]

excecutable = './TGF_Propa'

# loops over required initial parameters list
for _ in range(NB_copy_commands):
    for ini_alt in initial_alt_list:
        for beaming_ang in beaming_angle_list:
          for source_sigmat in source_sigma_time_list:
            for beaming_type in beaming_type_list:
               for tilt_ang in tilt_angle_list:
                   commands.append(excecutable +' '+str(nr_record_to_get_per_run)+' '+str(ini_alt)
                                   +' '+str(beaming_ang)+' '+ str(tilt_ang) +' '+str(beaming_type)+' '+str(source_sigmat))


# if number of commands is less than NB_CPU_LOCAL, fill the list with extra commands that are duplicate of previous ones
jj=0
while len(commands)<nprocs:
    commands.append(commands[jj])
    jj += 1
    if jj==len(commands):
        jj = 0

# print(commands)
command_number=len(commands)
print('Number of commands required '+ str(command_number))


####################################
computer_name = platform.node()

if "iftrom" in computer_name: # run on local (personal) computer
  
    nb_thread = NB_CPU_LOCAL # number of threads (cpu) to run
    
    # Making an array where each element is the list of command for a given thread

    pool = Pool(nb_thread) # to be always set to 1 for this MPI case
    for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
        if returncode != 0:
            print("%d command failed: %d" % (i, returncode))

else: # run on computer cluster using MPI

    listy = [[] for i in range(0,nprocs,1)]    
    
    i=0
    for j in range(0,command_number,1):
        listy[i].append(commands[j])
        i=i+1
        if i==nprocs:
            i=0

    pool = Pool(1) # to be always set to 1 for this MPI case
    for i, returncode in enumerate(pool.imap(partial(call, shell=True), listy[rank])):
        if returncode != 0:
            print("%d command failed: %d" % (i, returncode))

