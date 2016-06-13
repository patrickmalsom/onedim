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

parser.add_argument('--num', type=int, default=101,
         help='steps to perform in the simulation')
parser.add_argument('--dt', type=float, default=0.005,
         help='time step')
parser.add_argument('--eps', type=float, default=0.25,
         help='configurational temp')
parser.add_argument('--xstart', type=float, default=-0.4,
         help='starting position')
parser.add_argument('--method', type=str, default='midpt',
         help='quadrature: leapfrog, midpt, simpson')
         # leapfrog:0, midpt:1, simpson:2
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
clib.create_trajectory.restype = DOUBLE
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
  _fields_ = [('grn',DOUBLE)]
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
  

#class particle(xstart=-0.4):
#  def next_leapfrog(): 
traj=randGauss()

print(c_Pot(ctypes.c_double(.14)))

for i in range(params.num):
  traj[i].grn=np.random.rand()


print(c_create_trajectory(traj,params))
