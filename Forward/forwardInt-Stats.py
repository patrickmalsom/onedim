#!/bin/env python

# import libraries
import math
import sys
import time

import random
import numpy
numpy.random.seed(random.SystemRandom().randint(1,1000000))


##===================================
#import ctypes
## import the c code library
#clib=ctypes.CDLL("FatSkinny.so")
## library functions
#c_Force=clib.Force
#clib.Force.restype = ctypes.c_double
#clib.Force.argtypes = [ctypes.c_double]
#xold=(ctypes.c_double)()


##===================================
## set up the profiler
#import cProfile, pstats, StringIO
##profiler start
#pr = cProfile.Profile()
#pr.enable()

#===================================
# Define the constants
method=str(sys.argv[1])   # LeapFrog,MidPt,Simpson
Temp=float(sys.argv[2])   # configurational temperature
dt=float(sys.argv[3])     # time step dt
max_time=int(sys.argv[4]) # maximum time to run code (in seconds)

xstart=4/3.               # starting position
basin=xstart
LeftBasin=-2/3.
RightBasin=4/3.

#===================================
# define the force ( for this example V(x)=fat skinny )
def F(x): return x*(6.75 + x*(-5.0625 + x*(-11.390625 + (14.23828125 - 4.271484375*x)*x)))

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
    xguess = xold + dt* (0.16666666666666666*F(xold) + 0.16666666666666666*F(xnew) + 0.6666666666666666*F(0.5*(xold+xnew)) ) + noise
    while xnew-xguess>10.**-13.:
      xnew=xguess
      xguess = xold + dt* (0.16666666666666666*F(xold) + 0.16666666666666666*F(xnew) + 0.6666666666666666*F(0.5*(xold+xnew)) ) + noise
    return xguess

else:
  print "no method matching: " + method
  sys.exit(0)
#===================================

xold = xstart
pref=math.sqrt(2.0*Temp*dt)
signCt=math.copysign(1.,xstart)
trans=0
loop=0

start_time = time.time()

MaxLoops=sys.maxint-1

#===================================
while ( (time.time() - start_time) < max_time ) and (trans<10000) and (loop<MaxLoops):
  if abs(xold)>100. or math.isnan(xold):
    print "OVERFLOW: %s\t%f\t%f\t%d"%(method,Temp,dt,loop)
    sys.exit(0)
  for dummyLoop in xrange(10000):
    noise = pref*numpy.random.normal(0.0,1.0)

    xnew=genStep(xold,dt,noise)

    signCt+=math.copysign(1.0,xnew)
    if xnew>RightBasin and basin == RightBasin:
      basin = LeftBasin
    if xnew<LeftBasin and basin == LeftBasin:
      basin=RightBasin
      trans+=1

    xold=xnew
    loop+=1

##===================================
## print profile results
#pr.disable()
#s = StringIO.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
#print s.getvalue()

#===================================
print "%s\t%d\t%f\t%f\t%f\t%d"%(method,trans,Temp,dt,0.5*((signCt/loop)+1.0),loop)
