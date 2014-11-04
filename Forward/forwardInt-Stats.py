#!/usr/bin/env python

# import libraries
import math
import sys
import time

import random
import numpy
numpy.random.seed(random.SystemRandom().randint(1,1000000))

##===================================
## set up the profiler
#import cProfile, pstats, StringIO
##profiler start
#pr = cProfile.Profile()
#pr.enable()

##===================================
# Argparse: command line options
import argparse
parser = argparse.ArgumentParser(description='Forward integration of the SDE with MH(G)MC using Leap Frog, Mid Point and Simpsons.')
parser.add_argument('method', type=str, help='Integration method: LeapFrog MidPt Simpson')
parser.add_argument('MHMC', type=int, help='Metropolis-Hastings Boolean')
parser.add_argument('MHG', type=int, help='Metropolis-Hastings-Green Boolean')
parser.add_argument('eps', type=float, help='Configurational temperature')
parser.add_argument('dt', type=float, help='Time step')
parser.add_argument('NumB', type=int, help='Total path steps (number of beads)')
args = parser.parse_args()

method=args.method  # Integration method: LeapFrog,MidPt,Simpson
MHMC_bool=args.MHMC # Metropolis-Hastings switch (0:off, 1:on)
MHG_bool=args.MHG   # Metropolis-Hastings-Green switch (0:off, 1:on)
eps=args.eps        # Configurational temperature (eps)
dt=args.dt          # Time step (dt)
NumB=args.NumB      # Total path steps (number of beads)


##===================================
# save the decimal form of the method
if method == 'LeapFrog':
  method_decimal = 1
elif method == 'MidPt':
  method_decimal = 2
elif method == 'Simpson':
  method_decimal = 3
else:
  print "no method matching: " + method
  sys.exit(0)

#===================================
# define the potential, force and hessian ( for this example V(x)=fat skinny )
def U(x): return 1. + x*x*(-3.375 + x*(1.6875 + x*(2.84765625 + (-2.84765625 + 0.7119140625*x)*x)))
def F(x): return x*(6.75 + x*(-5.0625 + x*(-11.390625 + (14.23828125 - 4.271484375*x)*x)))
def H(x): return -6.75 + x*(10.125 + x*(34.171875 + x*(-56.953125 + 21.357421875*x)))


#===================================
# genStep: function to generate the xnew position
#    uses 'method' argument to determine integration technique
def genStep(int_method, xold, noise):
  # ============= Leap Frog Integrator ==============
  if int_method == "LeapFrog":
    return ( xold + dt*F(xold) + noise )

  # ============= Mid Point Integrator ==============
  elif int_method == "MidPt":

    # find the new position using leap frog
    xnew = xold + dt*F(xold) + noise

    # generate the new guess state using the mid point
    xguess = xold + dt*F(0.5*(xold+xnew)) + noise

    # iterate on this generation until the new step matches the guessed step to 13 digits
    while xnew-xguess>10.**-13.:
      xnew=xguess
      xguess = xold + dt*F(0.5*(xold+xnew)) + noise

    # perform 3 more of these guesses for good measure. This is important as the MHMC 
    # step is very sensitive to the correct future step being generated.
    for i in xrange(3):
      xnew=xguess
      xguess = xold + dt*F(0.5*(xold+xnew)) + noise
    
    # finally return the resulting future state
    return xguess

  # ============= Simpson's Integrator ==============
  elif int_method == "Simpson":
    # define dt/6.0
    dtOverSix = (dt/6.)

    # find the new position using leap frog
    xnew = xold + dt*F(xold) + noise

    # generate the new guess state using Simpson's method
    xguess = xold + dtOverSix*(F(xold) + F(xnew) + 4.*F(0.5*(xold+xnew)) ) + noise

    # iterate on this generation until the new step matches the guessed step to 13 digits
    while xnew-xguess>10.**-13.:
      xnew=xguess
      xguess = xold + dtOverSix*(F(xold) + F(xnew) + 4.*F(0.5*(xold+xnew)) ) + noise

    # perform 3 more of these guesses for good measure. This is important as the MHMC 
    # step is very sensitive to the correct future step being generated.
    for i in xrange(3):
      xnew=xguess
      xguess = xold + dtOverSix*(F(xold) + F(xnew) + 4.*F(0.5*(xold+xnew)) ) + noise

    # finally return the resulting future state
    return xguess

#===================================
# calcEnergy: function that calculates the energy error associated
#    with the integration forward. 
def calcEnergy(int_method, xold, xnew):
  # ============== Leap Frog Error =================
  if int_method == "LeapFrog":
    # evaluate the force at current and next step for later use
    Fx1, Fx0 = F(xnew), F(xold)
    # return the energy error. see notes for explanation of algebra.
    return 0.5*(xnew-xold)*(Fx1+Fx0)+dt*0.25*(Fx1*Fx1-Fx0*Fx0) + (U(xnew) - U(xold))

  # ============== Mid Point Error =================
  if int_method == "MidPt":
    # return the energy error. see notes for explanation of algebra.
    return (xnew-xold) * F(0.5*(xnew+xold)) + (U(xnew) - U(xold))

  # ============== Simpson's Error =================
  if int_method == "Simpson":
    # return the energy error. see notes for explanation of algebra.
    return 1/6. * (xnew-xold) * ( F(xnew)+F(xold)+4.*F(0.5*(xnew+xold)) ) + (U(xnew) - U(xold))
  
#===================================
# JacobianCorrection: function that calculates the correction to the MHMC method
#    to turn it into metropolis hastings green method. 

def JacobianCorrection(int_method, xold, xnew):
  # ============= Leap Frog Jacobian ===============
  if int_method == "LeapFrog":
    return 1.0
  # ============= Mid Point Jacobian ===============
  elif int_method == "MidPt":
    return 1.0
        
  # ============= Simpson's Jacobian ===============
  elif int_method == "Simpson":
    dtOverSix = (dt/6.0)
    xbar=0.5*(xold+xnew)
    return ( 1.0 + dtOverSix*(H(xold) + 2.0*H(xbar)) ) / ( 1.0 + dtOverSix*(H(xnew) + 2.0*H(xbar)) )
    
#===================================
# createTrajectory: main function that performs the simulation.
#    prints results to stdout when finished
def createTrajectory():

  # starting position
  xstart=4/3.
  # initial basin position
  basin=xstart

  # set the starting position 
  xold = xstart

  # prefactor in front of the noise to make v0h
  pref=math.sqrt(2.0*eps*dt)
  
  # counter that incriments if you are in pos or neg basin
  # Note: only works for barrier centered about 0.
  signCt=0.0

  # accumulators
  trans=0
  acc=0
  rej=0

  # basin minima
  LeftBasin=-2/3.
  RightBasin=4/3.

  #============== MAIN LOOP =====================
  for i in xrange(NumB):
    # generate thermalized noise to use when generating the new step
    v0h = pref*numpy.random.normal(0.0,1.0)

    # generate the new step (using specific method)
    xnew=genStep(method, xold, v0h)

    # incriment if the particle is in the right or left well
    signCt+=math.copysign(1.0,xnew)

    # counter for left/right basin. used to find total transition number
    if xnew>RightBasin and basin == RightBasin:
      basin = LeftBasin # dont incriment the transition counter
    if xnew<LeftBasin and basin == LeftBasin:
      basin=RightBasin
      trans+=1 # incriment the trasition counter. trans is back and forth. 

    # If metropolis hastings is required do the following loop
    if MHMC_bool==True:
      # ============= Metropolis Hastings Green =================
      # if metropolis-hastings-green is to be performed (includes the jacobian)
      if MHG_bool == True:
        if JacobianCorrection(method,xold,xnew)*math.exp(-calcEnergy(method,xold,xnew)/eps)>numpy.random.random():
          acc+=1
          xold=xnew
        else:
          rej+=1    
      # ============= Metropolis Hastings =================
      else:
        if math.exp(-calcEnergy(method,xold,xnew)/eps)>numpy.random.random():
          acc+=1
          xold=xnew
        else:
          rej+=1    

    # ============= NO Metropolis Test ===============
    else:
      xold=xnew


  # string to print at the end of the numeric part. Gives a human readable string
  endingPrintStr=method
  if MHMC_bool==True:
    endingPrintStr = endingPrintStr + "MH"
  if MHG_bool==True:
    endingPrintStr = endingPrintStr + "G"

  # print out the final results of the simulation
  print"%d %d %d %f %f %d %d %f %f #%s" % (method_decimal, MHMC_bool, MHG_bool, eps, dt, NumB, trans, 0.5*((signCt/NumB)+1.0), float(acc)/float(NumB), endingPrintStr)

# ============== run the code =======================
createTrajectory()
