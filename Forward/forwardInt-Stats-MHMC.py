#!/bin/env python

# import libraries
import math
import sys
import time

import random
import numpy
numpy.random.seed(random.SystemRandom().randint(1,1000000))

# Define the constants
method=sys.argv[1]   # LeapFrog,MidPt,Simpson
#Temp=0.15           # configurational temperature
#dt=0.025            # time step dt
#Nb=10000001           # path length


#===================================
# define the force ( for this example V(x)=fat skinny )
def F(x): return x*(6.75 + x*(-5.0625 + x*(-11.390625 + (14.23828125 - 4.271484375*x)*x)))

def U(x): return 1. + x*x*(-3.375 + x*(1.6875 + x*(2.84765625 + (-2.84765625 + 0.7119140625*x)*x)))

def H(x): return -6.75 + x*(10.125 + x*(34.171875 + x*(-56.953125 + 21.357421875*x)))
#===================================
# genStep: function to generate the xnew step under seperate methods
# Leap Frog integration
if method == "LeapFrog":
  def genStep(xold, dt, noise):
    return ( xold + dt*F(xold) + noise )

# Mid-Point integration
elif method == "MidPt":
  def genStep(xold, dt, noise):
    xnew = xold + dt*F(xold) + noise
    xguess = xold + dt*F(0.5*(xold+xnew)) + noise
    while xnew-xguess>10.**-13.:
      xnew=xguess
      xguess = xold + dt*F(0.5*(xold+xnew)) + noise
    return xguess

# Simpsons method integration
elif method == "Simpson":
  def genStep(xold, dt, noise):
    xnew = xold + dt*F(xold) + noise
    xguess = xold + (dt/6.)*(F(xold) + F(xnew) + 4.*F(0.5*(xold+xnew)) ) + noise
    while xnew-xguess>10.**-13.:
      xnew=xguess
      xguess = xold + (dt/6.)*(F(xold) + F(xnew) + 4.*F(0.5*(xold+xnew)) ) + noise
    for i in xrange(3):
      xnew=xguess
      xguess = xold + (dt/6.)*(F(xold) + F(xnew) + 4.*F(0.5*(xold+xnew)) ) + noise
    return xguess

else:
  print "no method matching: " + method
  sys.exit(0)
#===================================
#calcEnergy: function that calculates the energy error associated
#    with the integration forward. 
# Leap Frog Error
if method == "LeapFrog":
    def JacobianCorrection(xold,xnew,dt):
        return 1.0

    def calcEnergy(xold,xnew,dt):
        Fx1=F(xnew)
        Fx0=F(xold)
        return 0.5*(xnew-xold)*(Fx1+Fx0)+dt*0.25*(Fx1*Fx1-Fx0*Fx0)+U(xnew)-U(xold)

elif method == "MidPt":
    def JacobianCorrection(xold,xnew,dt):
        return 1.0

    def calcEnergy(xold,xnew,dt):
        return (xnew-xold)*F(0.5*(xnew+xold))+U(xnew)-U(xold)
        
elif method == "Simpson":
    def JacobianCorrection(xold,xnew,dt):
        oneSixth=0.16666666666666666
        oneThird=0.33333333333333333
        xbar=0.5*(xold+xnew)
        return (1. + dt*(oneSixth*H(xold) + oneThird*H(xbar))) / (1. + dt*(oneSixth*H(xnew) + oneThird*H(xbar)) )

    def calcEnergy(xold,xnew,dt):
        return 0.16666666666666666*(xnew-xold)*(F(xnew)+F(xold)+4.*F(0.5*(xnew+xold))) + (U(xnew) - U(xold))
    
    

    
def createTrajectory(temperature, dt, Nb, MHMC_bool):

    xstart=4/3.         # starting position
    basin=xstart
    LeftBasin=-2/3.
    RightBasin=4/3.
    
    xold = xstart
    pref=math.sqrt(2.0*temperature*dt)
    signCt=math.copysign(1.,xstart)
    trans=0
    acc=0
    rej=0

    #===================================
    for i in xrange(Nb):
        noise = pref*numpy.random.normal(0.0,1.0)

        xnew=genStep(xold,dt,noise)

        signCt+=math.copysign(1.0,xnew)
        if xnew>RightBasin and basin == RightBasin:
            basin = LeftBasin
        if xnew<LeftBasin and basin == LeftBasin:
            basin=RightBasin
            trans+=1

        if MHMC_bool==True:
            if JacobianCorrection(xold,xnew,dt)*math.exp(-calcEnergy(xold,xnew,dt)/temperature)>numpy.random.random():
                acc+=1
                xold=xnew
            else:
                rej+=1    
        else:
            xold=xnew

    print "JacobianNew %s  %d  %d  %f  %f  %f  %d  %f"%(method, MHMC_bool, trans,temperature,dt,0.5*((signCt/Nb)+1.0),Nb,float(acc)/float(Nb))

createTrajectory(0.15,0.025,1000000001,int(sys.argv[2]))
