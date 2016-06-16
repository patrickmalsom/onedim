#!/usr/bin/env python3

import numpy as np
import math
import sys
import os

# ===============================================
# ArgParse setup
# ===============================================
import argparse
parser = argparse.ArgumentParser(
        description='Forward B(s) statistics',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required Argument(s)
parser.add_argument('method', type=str,
         help='quadrature: leapfrog, midpt, simpson')
         # leapfrog:0, midpt:1, simpson:2

# Optional arguments
parser.add_argument('--num', type=int, default=101,
         help='steps to perform in the simulation')
parser.add_argument('--dt', type=float, default=0.005,
         help='time step, dt')
parser.add_argument('--eps', type=float, default=0.25,
         help='configurational temp')
parser.add_argument('--xstart', type=float, default=-0.4,
         help='starting position, x(0)')
parser.add_argument('--MHMC', type=int, default=0,
         help='Use MHMC boolean, 0:No 1:Yes')

#sys.argv includes a list of elements starting with the program
if len(sys.argv) < 2:
    parser.print_usage()
    sys.exit(0)

args = parser.parse_args()

# ===============================================
# Ctypes
# ===============================================
import ctypes
# define some c variable types
DOUBLE=ctypes.c_double
PDOUBLE=ctypes.POINTER(DOUBLE)
INT=ctypes.c_int
PINT=ctypes.POINTER(INT)

# import the c code library
if os.path.isfile("./calc_move.so") is False:
    print("ERROR: C file does not exist! Exiting...")
    sys.exit(1)
clib=ctypes.CDLL("./calc_move.so")

# Save the c functions to alternative names
c_Pot = clib.Pot
clib.Pot.restype = DOUBLE
c_create_trajectory = clib.create_trajectory
clib.create_trajectory.restype = INT
# returns acceptances for MHMC

# ===============================================
# Structs to pass to the C routines
# ===============================================
# NOTE: The fields MUST EXACTLY MATCH the structs in the C routine
# save parameters to a struct to pass to the C routines
class Parameters(ctypes.Structure):
  _fields_ = [('dt',DOUBLE),
              ('eps', DOUBLE),
              ('noisePref', DOUBLE),
              ('xstart', DOUBLE),
              ('num', INT),
              ('method', INT),
              ('MHMC', INT),
             ]
# Array struct
class traj_array(ctypes.Structure):
  _fields_ = [('grn',DOUBLE),
              ('pos',DOUBLE),
              ('rand',DOUBLE)
             ]
traj_array_type=traj_array*args.num

# ===============================================
# Initializing the parameters in the structs
# ===============================================
# Make instance of the Parameters class 
params=Parameters()
# Initialize the instance
params.dt = args.dt
params.num= args.num
params.eps= args.eps
params.xstart= args.xstart
params.noisePref=math.sqrt(2.0 * args.dt * args.eps)
if args.method == "leapfrog":
  params.method=0
elif args.method == "midpt":
  params.method=1
elif args.method == "simpson":
  params.method=2
else:
  sys.exit(0)
params.MHMC= args.MHMC
  
# ===============================================
# Function definitions
# ===============================================

def broad_frac(traj,params):
  # calculate the broad well fraction of the currently saved trajectory
  return sum( [ (np.sign(traj[i].pos)+1.0)*0.5 for i in  range(params.num)] )/params.num





# create an instance of the trajectory 
traj=traj_array_type()

# fill the Gaussian random numbers used as noise
for i in range(params.num):
  traj[i].grn=np.random.normal()

# if using MHMC, fill the random array used as acc/rej criteria in MHMC test
if params.MHMC == 1:
  for i in range(params.num):
    traj[i].rand = np.random.rand()


# generate a new trajectory
acceptances=c_create_trajectory(traj,params)
print(acceptances)

#for i in range(params.num):
#  print(traj[i].pos,traj[i].rand, traj[i].grn)
print(broad_frac(traj,params))
#print([ (np.sign(traj[i].pos)+1.0)*0.5 for i in  range(params.num)] )
