#!/bin/env python

# ====================================
# import libraries
import time
import math
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call

import ctypes

GenPathC=ctypes.CDLL("/home/patrick/git/onedim/guided/GenPaths.so")
CGenPaths=GenPathC.GenPaths


# ====================================
RunName="testa"

#m0
mPlus=0.969624
gammam=2.0
tm=5.

# A0
APlus=7.52137
AHalf=-7.0
gammaA=10.0
tA=5.0

# define some constants of the simulation
Nb=10001
Nu=Nb-1
NumPaths=10001
eps=0.15
xStart=0.969624
deltat=0.001
h=math.sqrt(2.0*deltat)
T=Nu*deltat
pref=math.sqrt(2.0*eps*deltat)

# ====================================
def mZero(t):
  return mPlus*math.tanh(gammam * (t-tm))

def AZero(t):
  return APlus + AHalf*math.exp(-gammaA*(t-tA)**2)

#def GenPaths(numPaths,mlist,Alist,RNGSeed):
#  mSaveName="mlist.dat"
#  ASaveName="Alist.dat"
#  np.savetxt(mSaveName,mlist)
#  np.savetxt(ASaveName,Alist)
#  call(["./guidedGenPaths.out",str(numPaths),str(RNGSeed),RunName,mSaveName,ASaveName])
#  return np.loadtxt("outfiles/"+RunName+"-"+str(numPaths)+".dat")
# ====================================
# make the control function lists
mlist=[ mZero(t) for t in np.linspace(0,10,Nb)]
Alist=[ AZero(t) for t in np.linspace(0,10,Nb)]

# ====================================

DOUBLE=ctypes.c_double
PDOUBLE=ctypes.POINTER(DOUBLE)
INT=ctypes.c_int


class Cbead(ctypes.Structure):
  _fields_ = [('m', DOUBLE),
              ('A', DOUBLE),
              ('xbar', DOUBLE),
              ('xxbar', DOUBLE),
              ('expVal', DOUBLE*5)
             ]

bead=Cbead*Nb

CGenPaths.argtypes = [INT,bead]
CGenPaths.restype = None

bead_Array = bead()

for i in range(Nb):
  bead_Array[i].m=mlist[i]
  bead_Array[i].A=Alist[i]
  bead_Array[i].xbar=0.0
  bead_Array[i].xxbar=0.0
  for j in range(5):
    bead_Array[i].expVal[j]=0.0


 
#doubleList= ctypes.c_double * Nb
#mlistC=doubleList(mlist)
#AlistC=doubleList(Alist)

CGenPaths(ctypes.c_int(10000),bead_Array)

plt.plot([bead_Array[i].xbar for i in range(Nb)])
plt.show()



#mlistC=np.array(mlist,dtype=ctypes.c_double)
#AlistC=np.array(Alist,dtype=ctypes.c_double)

#klstats=GenPaths(NumPaths,mlist,Alist,10)

# ====================================
# PLOTS
# ====================================

#timePlt=[t for t in np.linspace(0,T,Nb)]
#f, axarr = plt.subplots(2,2)
#
## A plot
#axarr[0,0].plot(timePlt, Alist)
#axarr[0,0].axis([0, 10, 0, 10])
#
## m plot
#axarr[0,1].plot(timePlt, mlist)
#axarr[0,1].plot(timePlt, 1.0/NumPaths*(klstats[:,2]))
#
## dD/dA plot
#axarr[1,0].plot(timePlt, 1.0/NumPaths*((klstats[:,4]*klstats[:,7]) - klstats[:,8]) )
#
## dD/dm plot
#axarr[1,1].plot(timePlt, 1.0/NumPaths*((klstats[:,4]*klstats[:,5]) - klstats[:,6]) )
#
#plt.show()
