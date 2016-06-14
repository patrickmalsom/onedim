#!/usr/bin/env python3

import numpy as np
import math
import sys
import os

import argparse
parser = argparse.ArgumentParser(
        description='Forward B(s) statistics',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
# Required Argument(s)
parser.add_argument('method', type=str,
         help='quadrature: leapfrog, midpt, simpson')
         # leapfrog:0, midpt:1, simpson:2

parser.add_argument('--num', type=int, default=101,
         help='steps to perform in the simulation')
parser.add_argument('--dt', type=float, default=0.005,
         help='time step')
parser.add_argument('--eps', type=float, default=0.25,
         help='configurational temp')
parser.add_argument('--xstart', type=float, default=-0.4,
         help='starting position')

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
if os.path.isfile("/home/patrick/forward_paths/calc_move.so") is False:
    print("ERROR: C file does not exist! Exiting...")
    sys.exit(1)
clib=ctypes.CDLL("/home/patrick/forward_paths/calc_move.so")

c_create_trajectory = clib.create_trajectory
c_Pot = clib.Pot
clib.Pot.restype = DOUBLE

# save parameters to a struct to pass to the C routines
class Parameters(ctypes.Structure):
  _fields_ = [('dt',DOUBLE),
              ('eps', DOUBLE),
              ('noisePref', DOUBLE),
              ('xstart', DOUBLE),
              ('num', INT),
              ('method', INT)
             ]

class traj_array(ctypes.Structure):
  _fields_ = [('grn',DOUBLE),
              ('pos',DOUBLE)
             ]
randGauss=traj_array*args.num

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
  
def broad_frac(traj,params):
  return sum( [ (np.sign(traj[i].pos)+1.0)*0.5 for i in  range(params.num)] )/params.num

#class particle(xstart=-0.4):
#  def next_leapfrog(): 
traj=randGauss()

for i in range(params.num):
  traj[i].grn=np.random.normal()

c_create_trajectory(traj,params)
#for i in range(params.num):
#  print(traj[i].pos)
print(broad_frac(traj,params))