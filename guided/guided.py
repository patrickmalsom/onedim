#!/bin/env python

# ==========================================================================
# Import Libraries
# ==========================================================================
import time
import math
import numpy as np
import matplotlib.pyplot as plt
import ctypes

# ==========================================================================
# Constants of the simulation
# ==========================================================================
Nb=10001
Nu=Nb-1
NumPaths=10001
eps=0.15
xStart=-0.969624
deltat=0.001
h=math.sqrt(2.0*deltat)
T=Nu*deltat
pref=math.sqrt(2.0*eps*deltat)

#m0
mPlus=0.969624
gammam=2.0
tm=5.

# A0
APlus=7.52137
AHalf=-7.0
gammaA=10.0
tA=5.0

# ==========================================================================
# C struct setup
# ==========================================================================
# import the c code library to generate paths
GenPathsLibrary=ctypes.CDLL("GenPaths.so")
# name the path generation function for convienence
CGenPaths=GenPathsLibrary.GenPaths

# define some c variable types
DOUBLE=ctypes.c_double
PDOUBLE=ctypes.POINTER(DOUBLE)
INT=ctypes.c_int
PINT=ctypes.POINTER(INT)

# make a class analogous to the struct in GenPaths.c
# this struct will be passed back with averages computed
class Averages(ctypes.Structure):
  _fields_ = [('m', DOUBLE),
              ('A', DOUBLE),
              ('xbar', DOUBLE),
              ('xxbar', DOUBLE),
              ('expVal', DOUBLE*5)
             ]

# define the array of averages that is to be passed to the C routine
beadType=Averages*Nb

class Parameters(ctypes.Structure):
  _fields_ = [('epsilon', DOUBLE),
              ('Numbead', INT),
              ('deltat', DOUBLE),
              ('xstart', DOUBLE)
             ]
# define the array of averages that is to be passed to the C routine
paramType=Parameters

# declare what arguments the function takes
CGenPaths.argtypes = [INT, INT, beadType]
CGenPaths.restype = None

# ==========================================================================
# Function declarations
# ==========================================================================
def mZero(t):
# base function for the mean: tanh
  return mPlus*math.tanh(gammam * (t-tm))

def AZero(t):
# base function for the width: gaussian well
  return APlus + AHalf*math.exp(-gammaA*(t-tA)**2)

def GenPaths(loops, RNGseed, mlist, Alist):
# call to the C library that will generate a path
# inputs:
#   loops: number of trajectories to calculate
#   mlist: input list for the controlling mean
#   Alist: input list for the controlling width
# outputs: 
#   returns an array with the expectation values (see wiki)

  # setting up the parameters
  params=paramType()
  params.epsilon=eps
  params.Numbead=Nb
  params.deltat=deltat
  params.xstart=xStart

  # allocate a new bead instance
  bead = beadType()
  # initialize the struct to pass to the c routine
  for i in range(Nb):
    bead[i].m=mlist[i]
    bead[i].A=Alist[i]
    bead[i].xbar=0.0
    bead[i].xxbar=0.0
    for j in range(5):
      bead[i].expVal[j]=0.0
  # call the C routine
  CGenPaths(ctypes.c_int(loops), ctypes.c_int(RNGseed), bead, params )
  # return a numpy array with the averages computed
  return np.array([ [bead[i].m,
                     bead[i].A,
                     bead[i].xbar,
                     bead[i].xxbar,
                     bead[i].expVal[0],
                     bead[i].expVal[1],
                     bead[i].expVal[2],
                     bead[i].expVal[3],
                     bead[i].expVal[4] 
                    ] for i in range(Nb) ])

# ==========================================================================
# Main routine
# ==========================================================================
# make the control function lists
mlist=[ mZero(t) for t in np.linspace(0,10,Nb)]
Alist=[ AZero(t) for t in np.linspace(0,10,Nb)]

#generate the paths and store stats in klstats array
klstats=GenPaths(NumPaths,13,mlist,Alist)

# ==========================================================================
# Plots!
# ==========================================================================
timePlt=[t for t in np.linspace(0,T,Nb)]
f, axarr = plt.subplots(2,2)

# A plot
axarr[0,0].plot(timePlt, Alist)
axarr[0,0].axis([0, 10, 0, 10])

# m plot
axarr[0,1].plot(timePlt, mlist)
axarr[0,1].plot(timePlt, 1.0/NumPaths*(klstats[:,2]))

# dD/dA plot
axarr[1,0].plot(timePlt, 1.0/NumPaths*((klstats[:,4]*klstats[:,7]) - klstats[:,8]) )

# dD/dm plot
axarr[1,1].plot(timePlt, 1.0/NumPaths*((klstats[:,4]*klstats[:,5]) - klstats[:,6]) )

plt.show()
