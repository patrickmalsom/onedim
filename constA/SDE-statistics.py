#!/bin/env python

import random
import math
from time import clock
import numpy as np
#import argparse
import sys

# =======================================================
# arguments for argparse error handling
# =======================================================
#input_args = argparse.ArgumentParser(description='Forward integration of SDE in an OU potential')
#input_args.add_argument('Method', metavar='Method', type=int, help='Integration method:    1:Euler 2:MidPt 3:Rotation')
#input_args.add_argument('MHMC', metavar='MHMC', type=int, help='perform MHMC:    1: YES   2: NO')
#input_args.add_argument('Potential', metavar='Potential', type=int, help='OU Potential:    1: 3->1.5   2: 10->1')
#input_args.add_argument('eps', metavar='eps', type=float, help='Configurational temperature (epsilon)')
#input_args.add_argument('xstart', metavar='xStart', type=float, help='Starting position')
#input_args.add_argument('deltat', metavar='deltat', type=float, help='Delta t')
#input_args.add_argument('Nb', metavar='Nb', type=int, help='Millions of Integration steps')
#args = input_args.parse_args()

# =======================================================
# Set vars with command line args
# =======================================================

# Method 1 -> Euler
# Method 2 -> Mid Point
# Method 3 -> Rotation
Method=int(sys.argv[1])

MHMC = int(sys.argv[2])

# Potential 1 -> A1=3 A2=1.5
# Potential 2 -> A1=10 A2=1
Potential=int(sys.argv[3])

eps=float(sys.argv[4])
xStart=float(sys.argv[5])
deltat = float(sys.argv[6])
Nb= int(sys.argv[7])


# =======================================================
# Define constants
# =======================================================

Nu=Nb-1
h=math.sqrt(2.0*deltat)
T= Nu * deltat
pref=math.sqrt(2.0*eps*deltat)

#make the list to store positions in
bins=np.array([0 for i in range(201)])
rejbins=np.array([0 for i in range(201)])


print("======================================")

if Method is 1:
  print("Euler     = %i" % Method)
if Method is 2:
  print("MidPt     = %i" % Method)
if Method is 3:
  print("Rotation  = %i" % Method)

if MHMC is 1:
  print("MHMC yes  = %i" % MHMC)
if MHMC is 2:
  print("MHMC no   = %i" % MHMC)

print("Potential = %i" % Potential)
print("Temp eps  = %f" % eps)
print("Start Pos = %f" % xStart)
print("Delta t   = %f" % deltat)
print("Num Steps = %i*10^6" % Nb)
print("########")


# =======================================================
# Potential 1
# =======================================================

if Potential is 1:
  def pot(x):
    if x < -0.10549438641447655:
      return 1.1060495262907208 + x*(2.6349336311983294 + 1.5692957403326768*x)
    elif x > 0.08435560927992086:
      return 1.0566826126123845 + (-1.8211246650661888 + 0.7846478701663384*x)*x
    else:
      return 1. + x*x*(-16.187357870924 + x*(12.336818433222824 + 324.3915486819529*x))

  def force(x):
    if x < -0.10549438641447655:
      return -2.6349336311983294 - 3.1385914806653537*x
    elif x > 0.08435560927992086:
      return 1.8211246650661888 - 1.5692957403326768*x
    else:
      return x*(32.374715741848 + (-37.01045529966847 - 1297.5661947278115*x)*x )

  def A(x):
    if x < -0.10549438641447655:
      return 3.1385914806653537
    elif x > 0.08435560927992086:
      return 1.5692957403326768
    else:
      return -32.374715741848 + x*(74.02091059933694 + 3892.698584183435*x)

# =======================================================
# Potential 2
# =======================================================

if Potential is 2:
  def pot(x):
    if x < -0.322174:
      return 1.7724051019087494 + x*(6.1648873109869555 + 5.360771574771266*x)
    elif x > 0.160554:
      return 1.0885716779019903 + (-1.5278198988098108 + 0.5360771574771266*x)*x
    else:
      return 0.9999999999999999 + x**2*(-7.395710880476022 + x*(6.600785187220504 + 30.72730601003735*x))

  def force(x):
    if x < -0.322174:
      return -6.1648873109869555 - 10.721543149542532*x
    elif x > 0.160554:
      return 1.5278198988098108-1.0721543149542532*x
    else:
      return x*(14.791421760952044 + (-19.802355561661514 - 122.9092240401494*x)*x)

  def A(x):
    if x < -0.322174:
      return 10.721543149542532
    elif x > 0.160554:
      return 1.0721543149542532
    else:
      return -14.791421760952044 + x*(39.60471112332303 + 368.7276721204482*x)

# =======================================================
# Potential 3
# =======================================================

if Potential is 3:
  def pot(x):
    return (3*x - 4)**4*(3*x + 2)**2/1024

  def force(x):
    return 27*x*(-81*x**4 + 270*x**3 - 216*x**2 - 96*x + 128)/512

  def A(x):
    return 10935*x**4/512 - 3645*x**3/64 + 2187*x**2/64 + 81*x/8 - 27/4

# =======================================================
# Sinc function
# =======================================================
def sinc(x2):
  if x2 > 0:
    sqx=math.sqrt(x2)
    return math.sin(sqx)/sqx
  elif x2<0:
    sqx=math.sqrt(-x2)
    return math.sinh(sqx)/sqx
  else:
    return 1.0

# =======================================================
# energy error
# =======================================================

if Method is 1:
  def deltaE(x0,x1):
    f0=force(x0)
    f1=force(x1)
    return pot(x1) - pot(x0) + 0.5*(x1-x0)*(f1+f0)+(h*h/8.0)*(f1*f1 - f0*f0)
    
if Method is 2:
  def deltaE(x0,x1):
    return pot(x1)-pot(x0)+(x1-x0)*force((x1+x0)/2.)
    
if Method is 3:
  def deltaE(x0,x1):
    return pot(x1)-pot(x0)+(x1-x0)*force((x1+x0)/2.)

# =======================================================
# Main Loop
# =======================================================
startTime = clock()

leftTotal=0.0
rightTotal=0.0
accTotal=0.0
rejTotal=0.0

trans=0
basin=np.sign(xStart)

x0=xStart
x1=xStart

for i in range(Nb):
  left=0
  right=0
  acc=0
  rej=0
  for j in range(1000000):
    # save the initial step and add it to the bin
    x0=x1
    bins[int((x0+10.)*10.)]+=1

    # generate the next step
    v0h=pref*random.gauss(0,1)
    x1=x0 + deltat*force(x0) + v0h
    
    if Method is 2:
      for i in range(4):
        guessx=x0+v0h+0.5*h*h*force((x0+x1)*0.5)
        x1=guessx
    
    if Method is 3:
      for i in range(4):
        phi2=2*A((x0+x1)*0.5)*deltat
        fourthphi2=phi2*0.25
        guessx=x0 + v0h *sinc(phi2) + deltat * sinc(fourthphi2)**2 * (force((x0+x1)*0.5) + A((x0+x1)*0.5) * (x1-x0)*0.5)
        x1=guessx

    if MHMC is 1:
      # do the  metropolis step
      if math.exp(-deltaE(x0,x1)/(eps))>random.random():
        acc=acc+1
      else:
        rej=rej+1
        x1=x0
        rejbins[int((x0+10.)*10.)]+=1
        
    else:
      acc=acc+1

    if abs(x1)>0.5 and x1*basin<0.0 :
      trans=trans+1
      basin=np.sign(x1)

    if x0 < 0:
      left=left+1
    else:
      right=right+1


  leftTotal=leftTotal+left/1000000.
  rightTotal=rightTotal+right/1000000.
  accTotal=accTotal+acc/1000000.
  rejTotal=rejTotal+rej/1000000.

elapsedTime = (clock() - startTime)

# =======================================================
# 
# =======================================================

print("Elapsed Time = %f" % elapsedTime)
print("left/right basin prob = %f" % (leftTotal/rightTotal) )
print("Left Total (*10-6) = %f" % (leftTotal) )
print("Right Total (*10-6) = %f" % (rightTotal) )
print("acc/rej = %f" % (accTotal/(accTotal+rejTotal)) )
print("Transitions = %i" % trans)
print("########")


print("Method, MHMC, Potential, eps, xStart, deltat, Nb, elapsedTime, leftTotal/rightTotal, leftTotal, rightTotal, (accTotal/(accTotal+rejTotal)), trans")
print("%i,%i,%i,%f,%f,%f,%i,%f,%f,%f,%f,%f,%i" %\
          (Method,MHMC,Potential,eps,xStart,deltat,Nb,elapsedTime,leftTotal/rightTotal,leftTotal,rightTotal,(accTotal/(accTotal+rejTotal)),trans) )

print np.array2string(np.linspace(-10,10,200), separator=', ',max_line_width=99999)[1:-1]
print np.array2string(bins, separator=', ',max_line_width=99999)[1:-1]
print np.array2string(rejbins, separator=', ',max_line_width=99999)[1:-1]
print("======================================")
