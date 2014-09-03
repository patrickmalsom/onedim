#!/usr/bin/env python2.7

# MALA routine

# ===============================================
# Import libraries
# ===============================================
import numpy as np
import math
import sys

# python profiler
import cProfile, pstats, StringIO

# ===============================================
# Ctypes
# ===============================================
import ctypes

# define some c variable types
DOUBLE=ctypes.c_double
PDOUBLE=ctypes.POINTER(DOUBLE)
INT=ctypes.c_int
PINT=ctypes.POINTER(INT)

# import the c code library
clib=ctypes.CDLL("NewMALA-func.so")

# library functions
#c_Pot=clib.Pot
#clib.Pot.restype = DOUBLE

c_calcForce=clib.calcForces
# fill out the Force variables for the struce

c_calcForcePrime=clib.calcForcePrimes
# fill out the ForcePrime variables for the struce

c_calcdg=clib.calcdg
# calculate dG/dx for the structure

c_calcSPDErhs=clib.calcSPDErhs
# calculate the right hand side vector for the SPDE

c_calcStateSPDE=clib.calcStateSPDE
# Function to calculate all of the structure parameters for the SPDE step

c_GaussElim=clib.GaussElim
# Gaussian elimination for computing L.x=b where L is the second deriv matrix

c_quadVar=clib.quadVar

# ===============================================
# constants/parameters
# ===============================================
deltat=0.005 # time step along the path
invdt=1/deltat # reciprocal of time step
eps=0.15 # configurational temperature
deltatau=0.00001 # time step between paths
noisePref=math.sqrt(4.0*eps*deltatau*invdt) # scaling for the random noise term
r=deltatau*0.5*invdt*invdt # constant matrix element in L (I+-rL)
NumB=100001 # number of 'beads' (total positions in path, path length)
NumU=10000 # number of time steps
NumL=99999 # size of matrix (I+-rL)
xStart=-2./3. # starting boundary condition
xEnd=4./3. # ending boundary condition

# ===============================================
# C struct setup
# ===============================================

# save parameters to a struct to pass to the C routines
class Parameters(ctypes.Structure):
  _fields_ = [('deltat',DOUBLE),
              ('invdt', DOUBLE),
              ('eps', DOUBLE),
              ('deltatau', DOUBLE),
              ('noisePref', DOUBLE),
              ('r', DOUBLE),
              ('NumB', INT),
              ('NumU', INT),
              ('NumL', INT),
              ('xStart', DOUBLE),
              ('xEnd', DOUBLE),
             ]
# define the array of averages that is to be passed to the C routine
paramType=Parameters

# make a class which is a clone of the C struct
# this structs pointer will be passed between C and Python code
class Averages(ctypes.Structure):
  _fields_ = [('pos', DOUBLE),
              ('randlist', DOUBLE),
              ('Force', DOUBLE),
              ('ForcePrime', DOUBLE),
              ('dg', DOUBLE),
              ('rhs', DOUBLE)
             ]
# define the array of averages that is to be passed to the C routine
pathType=Averages*NumB


# ===============================================
# Function declarations
# ===============================================
def Pot(x):
    return 1.+ x*x*(-3.375+x * (1.6875 +x * (2.84766 +(-2.84766+0.711914 * x) * x)))

def Force(x):
    return x * (6.75 + x * (-5.0625 + x * (-11.3906 + (14.2383 - 4.27148 * x) * x)))

def ForcePrime(x):
    return 6.75 + x * (-10.125 + x * (-34.1719 + (56.9531 - 21.3574 * x) * x))

def g(xm1,x0,x1):
    g1st=ForcePrime(x0)*Force((x0))
    g2nd=-0.5*invdt*(Force((x1))-2.0*Force((x0))+Force((xm1)))
    g3rd=-0.5*invdt*ForcePrime(x0)*(x1-2.0*x0+xm1)
    return g1st+g2nd+g3rd

def GaussElim(r,bVec):
    # Gaussian elimination for computing L.x=b where
    #   L is tridiagonal matrix
    #     mainDiag: 1+2r
    #     upper(lower)Diag: -r
    #   b is known vector: bVec is input
    #   x is unknown vector: returned at the end

    Num=len(bVec)
    
    al=[-r for i in range(Num-1)]
    am=[(1.0+2.0*r) for i in range(Num)]
    au=[-r for i in range(Num-1)]


    # Gaussian Elimination
    for i in range(Num-1):
        temp=-al[i]/am[i]
        #al[i]+=am[i]*temp
        am[i+1]+=au[i]*temp
        bVec[i+1]+=bVec[i]*temp

    # Back substitution
    for i in range(Num-1):
        temp=-au[Num-i-2]/am[Num-i-1]
        #au[Num-i-2]+=am[Num-i-1]*temp
        bVec[Num-i-2]+=bVec[Num-i-1]*temp
    

    # Divide by main diagonal
    for i in range(Num):
        bVec[i]=bVec[i]/am[i]
        #am[i]=am[i]/am[i]

    
    # Print the result
    return bVec


def rotatePaths():
    # rotate the structure pointers (old1->current0, current1->new0, new1->old0 )
    global pathOld,pathCur,pathNew
    temp=pathOld
    pathOld=pathCur
    pathCur=pathNew
    pathNew=temp

# ===============================================
# Unit Tests
def unitTest(path):
    # print the quadratic variation
    c_quadVar(path,params)
    # verify that the boundary conditions have not changed
    #if inpath[0] != outpath[0] and inpath[-1] != outpath[-1]:
    #    sys.exit("Boundary Conditions Changed! Aborting...")
    
def quadVar(poslist):
    quadVarSum=sum([(poslist[i]-poslist[i+1])**2 for i in range(len(poslist)-1)])
    quadVarSum/=(2.*eps*deltat*(len(poslist)-1))
    print "quad var: " + str(quadVarSum)

# check for the 
def checkBCs(path):
    return [path[0],path[-1]]

# ===============================================
# MAIN loop

np.random.seed(101)

inPath=np.loadtxt("SDEPath-fat-skinny.dat")

# ===============================================
# Initializing the structs
#declare the params struct
params=paramType()
params.deltat=deltat
params.invdt=invdt
params.eps=eps
params.deltatau=deltatau
params.noisePref=noisePref
params.r=r
params.NumB=NumB
params.NumU=NumU # not used
params.NumL=NumL

pathOld=pathType()
pathCur=pathType()
pathNew=pathType()


# testing temporary array
outPath=np.array([0.0 for i in np.arange(0,len(inPath),1)])

# start the profiler
pr = cProfile.Profile()
pr.enable()

for i in np.arange(0,len(inPath),1):
    pathOld[i].pos=inPath[i]


for loops in range(5):

    # noise vector
    for i in np.arange(1,len(inPath)-1,1):
        pathOld[i].randlist=np.random.normal(0,1)

    c_calcStateSPDE(pathOld,params)

    c_GaussElim(pathOld,pathCur,params)

    unitTest(pathCur)

    rotatePaths()


# profiler stop
pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print s.getvalue()

# test to see if the path is still correct
for i in np.arange(0,len(inPath),1):
    outPath[i]=pathOld[i].pos

if (np.loadtxt("outPathfive.dat") == outPath).all():
    print "SUCCESS!"
else:
    print "FAIL"

#np.savetxt("outPathC.dat",outPath)
