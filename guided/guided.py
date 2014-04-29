#!/usr/bin/env python2.7

# ==========================================================================
# Import Libraries
# ==========================================================================
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import time
import math
import ctypes
import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ==========================================================================
# Constants of the simulation
# ==========================================================================
# passed C params
Nb=10001
eps=0.15
xStart=-0.969624
deltat=0.001

# run specifics
NumPaths=1000
RNGseed=123

#m0
mPlus=0.969624
gammam=5.0
tm=5.0
# mH
gammamH=5.0
mHList=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# A0
APlus=7.52137
AHalf=-7.0
gammaA=10.0
tA=4.7
# AH
gammaAH=10.0
AHList=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#defined from above
Nu=Nb-1
h=math.sqrt(2.0*deltat)
T=Nu*deltat
pref=math.sqrt(2.0*eps*deltat)
#RNG stuff
rng.seed(RNGseed)

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

# ========== hermite polys ==========
def HH( t, n, tH, gamma):
    # normPref is 1/Sqrt[2^n n! Sqrt[Pi]]
    normPref = [0.751125544464942, 0.531125966013598, 0.265562983006799, 0.108415633823010, 0.0383307149314439, 0.0121212363525988]
    HHPoly = [ 1. , 0.  , 0.  , 0.   , 0. , 0. , \
               0. , 2.  , 0.  , 0.   , 0. , 0. , \
               -2., 0.  , 4.  , 0.   , 0. , 0. , \
               0. , -12., 0.  , 8.   , 0. , 0. , \
               12., 0.  , -48., 0.   , 16., 0. , \
               0. , 120., 0.  , -160., 0. , 32. ]
    tempSum = 0.0
    for k in range(6):
        tempSum +=  math.sqrt(gamma) * (t-tH)**k * HHPoly[(6*n) + k]
    return gamma**(0.25) * normPref[n] * tempSum

# ========== mean (m) ==========
def mZero(t):
# base function for the mean: tanh
  return mPlus*math.tanh(gammam * (t-tm))

def mHermite(t):
  return math.exp(-gammamH*(t-tm)**2) * np.sum( [ mHList[j] * HH(t,j,tm,gammamH)  for j in range(6) ])

def m(t):
  return mZero(t) + mHermite(t)

# ========== width (A) ==========
def AZero(t):
  return APlus + AHalf*math.exp(-gammaA*(t-tA)**2)

def AHermite(t):
  return math.exp(-gammaAH*(t-tA)**2) * np.sum( [ AHList[j] * HH(t,j,tA,gammaAH)  for j in range(6) ])

def A(t):
  return AZero(t) + AHermite(t)

# ========== calc paths ==========
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

# ========== approximate time lag ==========
def calcLagTime(mean, xbar, maxTimeLag):
  # input are the lists, m and xbar, and the maximum time shift necessary
  # output is an ordered pair list of the lag time tau and the d_L
  # $$ d_L=\int dt ( m(t) - x(t-\tau) )\^2 $$
  maxTauIndex=int(Nb*maxTimeLag/T)
  laglist=[ [0.0,0.0] for i in range(maxTauIndex)]
  print int((T - maxTimeLag)/deltat )
  for tau in range(maxTauIndex):
    dL=0.0
    for i in np.arange(maxTauIndex,Nb-maxTauIndex,1):
      dL+= (mean[i] - xbar[i+tau] )**2.0
    laglist[tau][0]=tau*deltat
    laglist[tau][1]=dL
  return laglist 

# ==========================================================================
# Main routine
# ==========================================================================

# generate m and A using the input parameters
mlist=np.array([m(t) for t in np.linspace(0,10,Nb)]);
Alist=np.array([A(t) for t in np.linspace(0,10,Nb)]);

for steepLoops in range(100):
  
  # generate the path stats using the C routine
  klstats=GenPaths(NumPaths,int(rng.randint(1,99999999)),mlist,Alist)
  normklstats=1.0/NumPaths
  #klstats
  #  0: mean
  #  1: A
  #  2: xbar
  #  3: xxbar
  #0    4: deltaU 
  #1    5: -A*(x-m)
  #2    6: -A*(x-m)*deltaU
  #3    7: 0.5*(x-m)^2
  #4    8: 0.5*(x-m)^2*deltaU
  templaglist=np.array(calcLagTime(mlist, normklstats*klstats[:,2], 0.4))
  #plt.plot(templaglist[:,0],templaglist[:,1])
  #plt.savefig("movies/HealAstrtest2.pdf", dpi = 100)
  #plt.close()

  dist=templaglist[:,1]
  lagderiv=[ (dist[i]-dist[i+1])/deltat for i in range(len(dist)-1)]

  # find the elemet where the deriv changes sign
  timeLag=np.where(np.diff(np.sign(lagderiv)))[0][0]
  print timeLag


  #plt.plot(templaglist[:,0][:-1],lagderiv)
  #plt.savefig("movies/HealAstrtest3.pdf", dpi = 100)
  #plt.close()

  polyFit= np.polyfit(np.array(templaglist[:,0]),np.array(templaglist[:,1]),2)
  print(polyFit)
  print(-2.0*polyFit[2]/polyFit[1])

  # do the mean gradient descent on the shifted kl deriv
  meanklstats = np.array( (normklstats*klstats[:,4]*normklstats*klstats[:,5])-normklstats*klstats[:,6] )
  print len(meanklstats)
  for i in range(Nb - timeLag):
    mlist[i] = mlist[i] - 0.00 * meanklstats[i+timeLag]

  # do the A gradient descent on the shifted kl deriv
  Aklstats = ( (normklstats*klstats[:,4]*normklstats*klstats[:,7])-normklstats*klstats[:,8] )
  for i in range(Nb - timeLag):
    Alist[i] = Alist[i] - 0.5 * Aklstats[i+timeLag]


  # ==========================================================================
  # Plots!
  # ==========================================================================
  timePlt=[t for t in np.linspace(0,T,Nb)]
  f, axarr = plt.subplots(2,2)

  # A plot
  axarr[0,0].plot(timePlt, Alist)
  #axarr[0,0].axis([0, 10, 0, 10])

  # m plot
  axarr[0,1].plot(timePlt, mlist)
  axarr[0,1].plot(timePlt, normklstats*(klstats[:,2]))
  #axarr[0,1].axis([4.9, 5.3, 0.0, 0.1])

  # dD/dA plot
  axarr[1,0].plot(timePlt, ((normklstats*klstats[:,4]*normklstats*klstats[:,7]) - normklstats*klstats[:,8]) )

  # dD/dm plot
  axarr[1,1].plot(timePlt, ((normklstats*klstats[:,4]*normklstats*klstats[:,5]) - normklstats*klstats[:,6]) )

  plt.savefig("movies/HealAstr"+str(100000+steepLoops)[1:]+".png", dpi = 100)
  plt.close(f)
#plt.show()
