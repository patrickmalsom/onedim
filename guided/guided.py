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

# import libraries
import random
import time
import math
import matplotlib.pyplot as plt
import numpy as np


# ================================================================
#                       Input Parameters
# ================================================================
# m0
mPlus=0.969624
gammam=2.0
tm=5.0
# A0
APlus=7.52137
AHalf=-7.0
gammaA=10.0
tA=5.0
# mH
gammamH=2.0
mHList=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# AH
gammaAH=10.0
AHList=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# define some constants of the simulation
Nb=10001
Nu=Nb-1
eps=0.15
xStart=-0.969624
deltat=0.001
h=math.sqrt(2.0*deltat)
T=Nu*deltat
pref=math.sqrt(2.0*eps*deltat)

# ================================================================
#                     Function Definitions
# ================================================================
def mZero(t):
    return mPlus*math.tanh(gammam * (t-tm))

def AZero(t):
    return APlus + AHalf*math.exp(-gammaA*(t-tA)**2)

def mHermite(t):
    
    # list of the nth hermite polynomial (max n is 15) and saves to a list 
    HermiteH = [\
    [1. , 0.  , 0.  , 0.   , 0. , 0. ],\
    [0. , 2.  , 0.  , 0.   , 0. , 0. ],\
    [-2., 0.  , 4.  , 0.   , 0. , 0. ],\
    [0. , -12., 0.  , 8.   , 0. , 0. ],\
    [12., 0.  , -48., 0.   , 16., 0. ],\
    [0. , 120., 0.  , -160., 0. , 32.]]


    # TODO need to figure out how to normalize these functions
    # computes the normalizations of the nth hermite polynomial (max n is 15) and saves to a list
    # HermiteHPref = [alpha**0.25 / math.sqrt(2.0**(n)*math.factorial(n)*math.sqrt(np.pi)) for n in range(16)]

    return math.exp(-gammamH*(t-tm)**2) * np.sum( [ mHList[j] * np.sum( [t**n * HermiteH[j][n]  for n in range(len(HermiteH))] ) for j in range(len(HermiteH)) ])

def m(t):
    return mZero(t) + mHermite(t)

def AHermite(t):
    
    # list of the nth hermite polynomial (max n is 15) and saves to a list 
    HermiteH = [\
    np.array([1. , 0.  , 0.  , 0.   , 0. , 0. ]),\
    np.array([0. , 2.  , 0.  , 0.   , 0. , 0. ]),\
    np.array([-2., 0.  , 4.  , 0.   , 0. , 0. ]),\
    np.array([0. , -12., 0.  , 8.   , 0. , 0. ]),\
    np.array([12., 0.  , -48., 0.   , 16., 0. ]),\
    np.array([0. , 120., 0.  , -160., 0. , 32.])]


    # TODO need to figure out how to normalize these functions
    # computes the normalizations of the nth hermite polynomial (max n is 15) and saves to a list
    # HermiteHPref = [alpha**0.25 / math.sqrt(2.0**(n)*math.factorial(n)*math.sqrt(np.pi)) for n in range(16)]

    return math.exp(-gammaAH*(t-tA)**2) * np.sum( [ AHList[j] * np.sum( [t**n * HermiteH[j][n]  for n in range(len(HermiteH))] ) for j in range(len(HermiteH)) ])

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

#initialize a few lists for storage of expectations
meanxPath=[0.0 for i in range(Nb)]
ev0=[0.0 for i in range(Nb)]
ev1=[0.0 for i in range(Nb)]
ev2=[0.0 for i in range(Nb)]
ev3=[0.0 for i in range(Nb)]
ev4=[0.0 for i in range(Nb)]

# number of loops per steepest descent loop
loops=1
# number of steep descent loops
SDloops=1


for runs in range(loops):
    xPath=[0.0 for i in range(Nb)]

    x0=xStart
    x1=xStart

    acc=0
    rej=0

    for i in range(Nb):
        x0=x1
        v0h=pref*random.gauss(0,1)
    
        t0=i*deltat
        t1=(i+1)*deltat
        thalf=(i+0.5)*deltat

        # TODO: is it mbar or m(tbar) for the mean???
        x1= genStep(x0, v0h, deltat, thalf)
        #x1= m0+((x0-m0)*(1.0 - 0.5*deltat*A0)+v0h)/(1+0.5*deltat*A0)
    
        deltaE = energyChange(x0, x1, t0, t1, thalf)

        if math.exp(-deltaE/eps) >= random.random():
            acc=acc+1

        else:
            rej=rej+1
            x1=x0
    
        xPath[i]=x1
        # calculate the expections and means
        meanxPath[i]+=x1
        mt1=m(t1)
        At1=A(t1)
        ev0[i]+= DeltaU(x1,mt1,At1)
        ev1[i]+= (-(x1-mt1)*At1)
        ev2[i]+= (-(x1-mt1)*At1)*DeltaU(x1,mt1,At1)
        ev3[i]+= (0.5*(x1-mt1)**2)
        ev4[i]+= (0.5*(x1-mt1)**2)*DeltaU(x1,mt1,At1)

# normalize the expectations
meanxPath= [ meanxPath[i]/float(loops) for i in range(Nb) ]
ev0= [ ev0[i]/float(loops) for i in range(Nb) ]
ev1= [ ev1[i]/float(loops) for i in range(Nb) ]
ev2= [ ev2[i]/float(loops) for i in range(Nb) ]
ev3= [ ev3[i]/float(loops) for i in range(Nb) ]
ev4= [ ev4[i]/float(loops) for i in range(Nb) ]


# =================================================================
#                     End of MAIN loop
# =================================================================

print time.time()-pretime

# <codecell>

print "accept: %d" % acc
print "reject: %d" % rej

#plt.plot([t for t in np.linspace(0,(Nb-1)*deltat,Nb)],xPath)
#plt.plot([t for t in np.linspace(0,(Nb-1)*deltat,Nb)],[m(t) for t in np.linspace(0,(Nb-1)*deltat,Nb)])
#plt.plot([t for t in np.linspace(0,(Nb-1)*deltat,Nb)],meanxPath)
#plt.plot([t for t in np.linspace(0,(Nb-1)*deltat,Nb)], [ 0.969624*math.tanh(2.0 * (t-5.)) for t in np.linspace(0,10,Nb)])
#plt.show()

