#!/usr/bin/env python

# Importing libraries
import random
import math
from time import clock
from scipy import integrate
import matplotlib.pyplot as plt

# Define constants
Nb= 100001 # maximum integration steps
Nu=Nb-1 #total time steps
eps=.25 # configurational temperature
deltat = .001 #step size
h=math.sqrt(2.0*deltat) # integration time step
T= Nu * deltat # total integration time
pref=math.sqrt(2.0*eps*deltat) # prefactor for thermal noise
xStart=-1. # simulation's starting position

# Define the potentials
def pot(x): 
  return (x**2.0 -1.0)**2.0
def force(x):
  return -4.0*x*(x**2.0 -1.)

savedPaths=0 # number of accepted sample paths
totalPaths=100
transOffset=2000 #integration steps to overshoot
avePath=[0.0 for i in range(transOffset*2)] # intialize an array to store the path
transMark=0.5 # position to count the transition at

while savedPaths<totalPaths:
  # Forward integration initialization
  x0=xStart # initialize x0: current step
  x1=xStart # initialize x1: future step
  switch=0 # counter: saves steps in the opposite well
  path=[0.0 for i in range(Nb)] # intialize an array to store the path
  index=0 # counter for current integration position (index of path array)
  
  # perform the forward integration
  # stop when the path changes basin
  while index<Nb and (x0<transMark or switch<transOffset):
    x0=x1
    if x0>transMark or switch>1:
      switch+=1
    v0h=pref*random.gauss(0,1)
    x1=x0 + deltat*force(x0) + v0h
    for i in range(4):
      guessx=x0+v0h+0.5*h*h*force((x0+x1)*0.5)
      x1=guessx
  
    path[index]=x1
    index+=1
  
  #find middle of trans trigger
  tempIndex = index
  posTrigger=tempIndex-transOffset
  while path[tempIndex-1] > -0.5:
    tempIndex-=1
  negTrigger=tempIndex
  middle=int(0.5*(posTrigger + negTrigger))

  print "index: %i i| pos: %f | saved: %i | %s %s %s " % (index, path[index-1], savedPaths, path[Nu] is 0.0, path[index-1]>0.0, index > transOffset*2)
  
  # decide to keep the path or not
  if path[Nu] is 0.0 and path[index-1]>0.0 and index > transOffset*2:
    # accepted path
    accPath=path[middle-transOffset:middle+transOffset]
    savedPaths+=1 #incriment the saved paths counter
    tempPath=avePath
    avePath=[accPath[i]+tempPath[i] for i in range(transOffset*2)]
  

tempPath=avePath
avePath=[tempPath[i]/float(totalPaths) for i in range(transOffset*2)]

plt.plot(avePath)
plt.show()
