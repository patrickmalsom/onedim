#!/usr/bin/env python

import numpy as np
import random
import math
import scipy
import scipy.integrate
import sys

# set the random seed of numpy
np.random.seed(random.SystemRandom().randint(1,1000000))

boltz= lambda x : scipy.exp(-4*((8 - 5*x)**8 * (2 + 5*x)**2)/(2**26))
Z0=scipy.integrate.quad(boltz,-np.infty,np.infty)[0]

def draw_from_pdf():
    
    # define functions
    eps=0.25
    potential = lambda x : ((8 - 5*x)**8 * (2 + 5*x)**2)/(2**26)
    boltz= lambda x : scipy.exp(-potential(x)/eps)
    Z0 = scipy.integrate.quad(boltz,-np.infty,np.infty)[0]
    PDF = lambda x : 1/Z0*math.exp(-potential(x)/eps)
    CDF = lambda a : scipy.integrate.quad(PDF,-np.infty,a)[0]
    
    guess=-5. # Value where the CDF is zero
    delta=1. # inital jump size
    
    random_number=np.random.rand()
    
    # calculate value of CDF at the guess value and test if rand<CDF(guess)
    while delta >0.000000001:
        while CDF(guess)<=random_number:
            guess+=delta
        guess-=delta
        delta/=2.
    return guess

pot = lambda y : (-5*y + 8)**8*(5*y + 2)**2/67108864
force = lambda y : -125*y*(5*y - 8)**7*(5*y + 2)/33554432

eps=0.25
dt=0.005
NumB=int(sys.argv[3])

# basin minima
LeftBasin=-2/5.
RightBasin=8/5.

# genStep: function to generate the xnew position
genStep = lambda xold, noise : ( xold + dt*force(xold) + noise )

def forwardTrajectory(xstart):
    traj=[0.0 for i in range(NumB)]

    # set the starting position 
    xold = xstart

    # prefactor in front of the noise to make v0h
    pref=math.sqrt(2.0*eps*dt)
    
    #============== MAIN LOOP =====================
    for i in range(NumB):
        # save position to trajectory list
        traj[i]=xold
        
        # generate thermalized noise to use when generating the new step 
        #   noise = v_0 * h = sqrt(2*dt) * sqrt(eps) * xi
        noise = pref*np.random.normal(0.0,1.0)

        
        # generate the new step (using specific method)
        xnew=genStep(xold, noise)
        
        xold=xnew
        
        
    return traj

# no pinning. start at drawn point and end wherever
if int(sys.argv[2]) == 0:
    for i in range(int(sys.argv[1])):
        print "nopin Nb%i %0.13f" % (  int(sys.argv[3]), np.sum((np.sign(forwardTrajectory(draw_from_pdf()))+1)/2)/NumB)

# start pinned at -2/5. end not pinned
if int(sys.argv[2]) == 1:
    for i in range(int(sys.argv[1])):
        print "startpin Nb%i %0.13f" % (  int(sys.argv[3]), np.sum((np.sign(forwardTrajectory(-2/5.))+1)/2)/NumB)

# start pinned at -2/5. end must be positive
if int(sys.argv[2]) == 2:
    i,j = 0,0
    while i < int(sys.argv[1]):
        temp_trajectory = forwardTrajectory(-2/5.)
        if temp_trajectory[-1]>0:
            print "bothpin Nb%i %0.13f" % ( int(sys.argv[3]),  np.sum((np.sign(temp_trajectory)+1)/2)/NumB)
            i+=1
        j+=1
    print "rejrate %i %i" %(i,j)

# start pinned at 8/5. end not pinned
if int(sys.argv[2]) == 3:
    for i in range(int(sys.argv[1])):
        print "startrevpin Nb%i %0.13f" % (  int(sys.argv[3]), np.sum((np.sign(forwardTrajectory(8/5.))+1)/2)/NumB)

# start pinned at 8/5. end must be negative
if int(sys.argv[2]) == 4:
    i,j = 0,0
    while i < int(sys.argv[1]):
        temp_trajectory = forwardTrajectory(8/5.)
        if temp_trajectory[-1]<0:
            print "bothrevpin Nb%i %0.13f" % ( int(sys.argv[3]),  np.sum((np.sign(temp_trajectory)+1)/2)/NumB)
            i+=1
        j+=1
    print "rejraterev %i %i" %(i,j)
