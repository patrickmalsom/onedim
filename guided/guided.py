#!/usr/bin/env python

# ==================================
# guided HMC in one dim
# Patrick Malsom, spring 2014
# ==================================

# import libraries
import random
import time
import math
import matplotlib.pyplot as plt
import numpy as np


# define some constants of the simulation
Nb=10001
Nu=Nb-1
eps=0.15
xStart=-0.969624
deltat=0.001
h=math.sqrt(2.0*deltat)
T=Nu*deltat
pref=math.sqrt(2.0*eps*deltat)

# list of the nth hermite polynomial (max n is 15) and saves to a list 
# (hardcoded for compatibility with older nupy versions)
#HermiteH = [np.array([1.]),\
#np.array([0.,2.]),\
#np.array([-2.,0.,4.]),\
#np.array([0.,-12.,0.,8.]),\
#np.array([12.,0.,-48.,0.,16.]),\
#np.array([0.,120.,0.,-160.,0.,32.]),\
#np.array([-120.,0.,720.,0.,-480.,0.,64.]),\
#np.array([0.,-1680.,0.,3360.,0.,-1344.,0.,128.]),\
#np.array([1680.,0.,-13440.,0.,13440.,0.,-3584.,0.,256.]),\
#np.array([0.,30240.,0.,-80640.,0.,48384.,0.,-9216.,0.,512.]),\
#np.array([-30240.,0.,302400.,0.,-403200.,0.,161280.,0.,-23040.,0.,1024.]),\
#np.array([0.00000000e+00,-6.65280000e+05,0.00000000e+00,2.21760000e+06,0.00000000e+00,-1.77408000e+06,0.00000000e+00,5.06880000e+05,0.00000000e+00,-5.63200000e+04,0.00000000e+00,2.04800000e+03]),\
#np.array([6.65280000e+05,0.00000000e+00,-7.98336000e+06,0.00000000e+00,1.33056000e+07,0.00000000e+00,-7.09632000e+06,0.00000000e+00,1.52064000e+06,0.00000000e+00,-1.35168000e+05,0.00000000e+00,4.09600000e+03]),\
#np.array([0.00000000e+00,1.72972800e+07,0.00000000e+00,-6.91891200e+07,0.00000000e+00,6.91891200e+07,0.00000000e+00,-2.63577600e+07,0.00000000e+00,4.39296000e+06,0.00000000e+00,-3.19488000e+05,0.00000000e+00,8.19200000e+03]),\
#np.array([-1.72972800e+07,0.00000000e+00,2.42161920e+08,0.00000000e+00,-4.84323840e+08,0.00000000e+00,3.22882560e+08,0.00000000e+00,-9.22521600e+07,0.00000000e+00,1.23002880e+07,0.00000000e+00,-7.45472000e+05,0.00000000e+00,1.63840000e+04]),\
#np.array([0.00000000e+00,-5.18918400e+08,0.00000000e+00,2.42161920e+09,0.00000000e+00,-2.90594304e+09,0.00000000e+00,1.38378240e+09,0.00000000e+00,-3.07507200e+08,0.00000000e+00,3.35462400e+07,0.00000000e+00,-1.72032000e+06,0.00000000e+00,3.27680000e+04])]
# computes the normalizations of the nth hermite polynomial (max n is 15) and saves to a list
#HermiteHPref = [alpha**0.25 / math.sqrt(2.0**(n)*math.factorial(n)*math.sqrt(np.pi)) for n in range(16)]

def m(t):
    return 0.969624*math.tanh(2.0 * (t-4.8))

def A( t):
    return -7.0*math.exp( -4.0 * (t - 5.0) * (t-5.0) )+7.52137

def U0(x, mean, width):
    return 0.5*( (x - mean) * (x - mean) ) * width

def U(x):
    return (x*x - 1.0)**2.0

def DeltaU(x, mean, width):
    return U(x) - U0(x,mean,width)

def genStep(x0, v0h, deltat, thalf):
    return (x0 + v0h - 0.5*deltat*x0*A(thalf) + m(thalf)*A(thalf)*deltat) / (1. + 0.5*deltat*A(thalf))

def energyChange(x0, x1, t0, t1, thalf):
    return U0(x1,m(t1),A(t1)) - U0(x0,m(t0),A(t0)) - (x1-x0)*A(thalf)*(x1+x0-2.0*m(thalf))*0.5

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
loops=10
# number of steep descent loops
SDloops=1


for runs in range(loops):
    xPath=[0.0 for i in range(Nb)]

    x0=xStart
    x1=xStart

    #m0=-0.969624
    #A0=7.52137

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
        ev0[i]+= DeltaU(x1,m(t1),A(t1))
        ev1[i]+= (-(x1-m(t1))*A(t1))
        ev2[i]+= (-(x1-m(t1))*A(t1))*DeltaU(x1,m(t1),A(t1))
        ev3[i]+= (0.5*(x1-m(t1))**2)
        ev4[i]+= (0.5*(x1-m(t1))**2)*DeltaU(x1,m(t1),A(t1))

# normalize the expectations
meanxPath= [ meanxPath[i]/float(loops) for i in range(Nb) ]
ev0= [ ev0[i]/float(loops) for i in range(Nb) ]
ev1= [ ev1[i]/float(loops) for i in range(Nb) ]
ev2= [ ev2[i]/float(loops) for i in range(Nb) ]
ev3= [ ev3[i]/float(loops) for i in range(Nb) ]
ev4= [ ev4[i]/float(loops) for i in range(Nb) ]

print time.time()-pretime

# <codecell>

print "accept: %d" % acc
print "reject: %d" % rej

#plot([t for t in linspace(0,(Nb-1)*deltat,Nb)],xPath)
#plot([t for t in linspace(0,(Nb-1)*deltat,Nb)],[m(t) for t in linspace(0,(Nb-1)*deltat,Nb)])
#plot([t for t in linspace(0,(Nb-1)*deltat,Nb)],meanxPath)
#plot([t for t in linspace(0,(Nb-1)*deltat,Nb)], [ 0.969624*math.tanh(2.0 * (t-5.)) for t in linspace(0,10,Nb)])
#show()

