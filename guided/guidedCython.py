#!/usr/bin/env python

# ==================================
# guided HMC in one dim
# Patrick Malsom, spring 2014
# ==================================

# Program outline
# input control functions: m(t) and A(t)
# input potential: U(t)
# input hermite params for m and A
#
# run the simulation many times and generate the expectation values ev{0,4}
# take the dot product dD/dA * dA/da * HH * {A,m}
# make changes to the hermite params using steepest descent
# repeat entire routine

# Import cython functions. Compile with:
#     python guidedCythonSetup.py build_ext --inplace
from guidedCython import *


# import libraries
import time
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ================================================================
#                       Class Definitions
# ================================================================
class params:
  """ class to store parameters"""
  mPlus = None
  gammam = None
  tm = None
  gammamH = None
  mHList = None
  APlus  = None
  AHalf = None
  gammaA = None
  tA = None
  gammaAH = None
  AHList = None

class stats:
  ev0 = None
  ev1 = None
  ev2 = None
  ev3 = None
  ev4 = None
  xpath = None
  xxpath = None
  
  
# ================================================================
#                       Input Parameters
# ================================================================
# m0
params.mPlus=0.969624
params.gammam=2.0
params.tm=5.
# mH
params.gammamH=2.0
params.mHList=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# for plots
mPrintString= "mean parameters \n m+: " + str(params.mPlus) + "\n gammam: " + str(params.gammam) + "\n tm: " + str(params.tm) + "\n gammamH: " + str(params.gammamH) + "\n " + str(params.mHList)

# A0
params.APlus=7.52137
params.AHalf=-7.0
params.gammaA=10.0
params.tA=5.0
# AH
params.gammaAH=10.0
params.AHList=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# for plots
APrintString= "A parameters \n A+: " + str(params.APlus) + "\n AHalf: " + str(params.AHalf) + "\n gammaA: " + str(params.gammaA) + "\n tA: " + str(params.tA) + "\n gammaAH: " + str(params.gammaAH) + "\n " + str(params.AHList)

# define some constants of the simulation
Nb=10001
Nu=Nb-1
eps=0.15
xStart=0.969624
deltat=0.001
h=math.sqrt(2.0*deltat)
T=Nu*deltat
pref=math.sqrt(2.0*eps*deltat)

# ================================================================
#                     Function Definitions
# ================================================================
#def HH(t,n,tH):
#    #returns the value of the nth hermite polynomial centered about time tH
#    HermitePolys = [\
#    np.array([1. , 0.  , 0.  , 0.   , 0. , 0. ]),\
#    np.array([0. , 2.  , 0.  , 0.   , 0. , 0. ]),\
#    np.array([-2., 0.  , 4.  , 0.   , 0. , 0. ]),\
#    np.array([0. , -12., 0.  , 8.   , 0. , 0. ]),\
#    np.array([12., 0.  , -48., 0.   , 16., 0. ]),\
#    np.array([0. , 120., 0.  , -160., 0. , 32.])]
#    # TODO need to figure out how to normalize these functions
#    # computes the normalizations of the nth hermite polynomial (max n is 15) and saves to a list
#    # HermiteHPref = [alpha**0.25 / math.sqrt(2.0**(n)*math.factorial(n)*math.sqrt(np.pi)) for n in range(16)]
#    return np.sum( [ (t-tH)**k * HermitePolys[n][k] for k in range(len(HermitePolys[n])) ] )

def mZero(t):
  return params.mPlus
#  return params.mPlus*math.tanh(params.gammam * (t-params.tm))

def mHermite(t):
  return math.exp(-params.gammamH*(t-params.tm)**2) * np.sum( [ params.mHList[j] * HH(t,j,params.tm,params.gammamH)  for j in range(6) ])

def m(t):
  return mZero(t) + mHermite(t)
#  return mZero(t) + mHermite(t)

def AZero(t):
  return params.APlus + params.AHalf*math.exp(-params.gammaA*(t-params.tA)**2)

def AHermite(t):
  return math.exp(-params.gammaAH*(t-params.tA)**2) * np.sum( [ params.AHList[j] * HH(t,j,params.tA,params.gammaAH)  for j in range(6) ])

def A(t):
  return AZero(t) + AHermite(t)

def U(x):
  return (x*x - 1.0)**2.0

# =================================================================
#                     Static Functions
# ================================================================
def U0(x, mean, width):
  return 0.5*( (x - mean) * (x - mean) ) * width

def DeltaU(x, mean, width):
  return U(x) - U0(x,mean,width)

def genStep(x0, v0h, deltat, thalf):
  return (x0 + v0h - 0.5*deltat*x0*A(thalf) + m(thalf)*A(thalf)*deltat) / (1. + 0.5*deltat*A(thalf))

def energyChange(x0, x1, t0, t1, thalf):
  return U0(x1,m(t1),A(t1)) - U0(x0,m(t0),A(t0)) - (x1-x0)*A(thalf)*(x1+x0-2.0*m(thalf))*0.5


# =================================================================
#                     Start of MAIN loop
# =================================================================

# Start the timer
pretime=time.time()


# number of loops per steepest descent loop
loops=10
# number of steep descent loops
SDloops=1

for descents in range(SDloops):

  #initialize a few lists for storage of expectations
  meanxPath=[0.0 for i in range(Nb)]
  stats.ev0=[0.0 for i in range(Nb)]
  stats.ev1=[0.0 for i in range(Nb)]
  stats.ev2=[0.0 for i in range(Nb)]
  stats.ev3=[0.0 for i in range(Nb)]
  stats.ev4=[0.0 for i in range(Nb)]

  for runs in range(loops):
    xPath=[0.0 for i in range(Nb)]

    x0=xStart
    x1=xStart

    acc=0
    rej=0

    for i in range(Nb):
      x0=x1
      v0h=pref*np.random.normal(0,1)

      t0=i*deltat
      t1=(i+1)*deltat
      thalf=(i+0.5)*deltat

      # TODO: is it mbar or m(tbar) for the mean???
      x1= genStep(x0, v0h, deltat, thalf)
      #x1= m0+((x0-m0)*(1.0 - 0.5*deltat*A0)+v0h)/(1+0.5*deltat*A0)

      deltaE = energyChange(x0, x1, t0, t1, thalf)

      if math.exp(-deltaE/eps) >= np.random.random():
        acc=acc+1

      else:
        rej=rej+1
        x1=x0

      xPath[i]=x1
      # calculate the expections and means
      meanxPath[i]+=x1
      mt1=m(t1)
      At1=A(t1)
      stats.ev0[i]+= DeltaU(x1,mt1,At1)
      stats.ev1[i]+= (-(x1-mt1)*At1)
      stats.ev2[i]+= (-(x1-mt1)*At1)*DeltaU(x1,mt1,At1)
      stats.ev3[i]+= (0.5*(x1-mt1)**2)
      stats.ev4[i]+= (0.5*(x1-mt1)**2)*DeltaU(x1,mt1,At1)

  # normalize the expectations
  meanxPath= np.array([ meanxPath[i]/float(loops) for i in range(Nb) ] )
  stats.ev0= np.array([ stats.ev0[i]/float(loops) for i in range(Nb) ] )
  stats.ev1= np.array([ stats.ev1[i]/float(loops) for i in range(Nb) ] )
  stats.ev2= np.array([ stats.ev2[i]/float(loops) for i in range(Nb) ] )
  stats.ev3= np.array([ stats.ev3[i]/float(loops) for i in range(Nb) ] )
  stats.ev4= np.array([ stats.ev4[i]/float(loops) for i in range(Nb) ] )


  # =================================================================
  #                     Gradient Descent
  # =================================================================
  # This is more complex than normal as there are many different parameters to minimize.
  # hermite part of A
  print "mean gradients:"
  meanGrads=np.array( [ np.sum( ( (stats.ev0*stats.ev1)-stats.ev2 ) * np.array([ math.exp(-params.gammamH *(t-params.tm)**2.0) * HH(t,n,params.tm,params.gammamH) for t in np.linspace(0,10,Nb) ]) ) for n in range(6)] )
  print list(meanGrads)
  print "new mean parameters:"
  params.mHList=list(np.array( params.mHList ) + 0.0*np.array([0.1,0.1,0.1,0.1,0.1,0.1]) * meanGrads)
  print params.mHList
  print " "
  
  print "A gradients:"
  AGrads=np.array( [ np.sum( ( (stats.ev0*stats.ev3)-stats.ev4 ) * np.array([ math.exp(-params.gammamH *(t-params.tA)**2.0) * HH(t,n,params.tm,params.gammaAH) for t in np.linspace(0,10,Nb) ]) ) for n in range(6)] )
  print list(AGrads)
  print "new A parameters:"
  params.AHList = list(np.array( params.AHList ) - 0.005*np.array([1.0,1.0,1.0,1.0,1.0,1.0]) * AGrads)
  print params.AHList
  print " "
  
  
  print time.time()-pretime
  
  # <codecell>
  
  print "accept: %d" % acc
  print "reject: %d" % rej
  
  # save the A as a picture
#  plt.plot([t for t in np.linspace(0,T,Nb)], [A(t) for t in np.linspace(0,T,Nb)])
#  plt.axis([0, 10, 0, 10])


  # =================================================================
  #                         Plots!
  # =================================================================
  
  timePlt=[t for t in np.linspace(0,T,Nb)]
  f, axarr = plt.subplots(2,2)

  # A plot
  axarr[0,0].plot(timePlt, [A(t) for t in np.linspace(0,T,Nb)])
  axarr[0,0].axis([0, 10, 0, 10])
  #axarr[0,0].text(6., 1., APrintString)
  
  # m plot
  axarr[0,1].plot(timePlt, [m(t) for t in np.linspace(0,(Nb-1)*deltat,Nb)])
  axarr[0,1].plot(timePlt, meanxPath)
  #axarr[0,1].text(6., -1., mPrintString)
  #axarr[0,1].plot(timePlt, [ 0.969624*math.tanh(2.0 * (t-5.)) for t in np.linspace(0,10,Nb)])
  
  # dD/dA plot
  axarr[1,0].plot(timePlt, (stats.ev0*stats.ev3) - stats.ev4 )

  # dD/dm plot
  axarr[1,1].plot(timePlt, (stats.ev0*stats.ev1) - stats.ev2 )
  
  plt.savefig("HealAstr"+str(1000+descents)[1:]+".jpg", dpi = 400)
#  plt.show()
  

