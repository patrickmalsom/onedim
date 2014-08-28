#!/usr/bin/env python2.7

# MALA routine

# ===============================================
# Import libraries
import numpy as np
import math
import sys
import time

# ==========================================================================
# Ctypes
# ==========================================================================
import ctypes
from ctypes import c_double
from ctypes import c_int

# define some c variable types
DOUBLE=ctypes.c_double
PDOUBLE=ctypes.POINTER(DOUBLE)
INT=ctypes.c_int
PINT=ctypes.POINTER(INT)
# import the c code library
clib=ctypes.CDLL("NewMALA-func.so")
# library functions (named for cinvienence)
#c_Pot=clib.Pot
#clib.Pot.restype = ctypes.c_double

c_calcForce=clib.calcForce
c_calcForcePrime=clib.calcForcePrime
c_calcdg=clib.calcdg

c_calcSPDErhs=clib.calcSPDErhs
c_calcStateSPDE=clib.calcStateSPDE

c_GaussElim=clib.GaussElim


# ===============================================
# constants/parameters
# ===============================================
deltat=0.005
invdt=1/deltat
eps=0.15
deltatau=0.00001
noisePref=math.sqrt(4.0*eps*deltatau*invdt)
r=deltatau*0.5*invdt*invdt
NumB=100001
NumU=10000
NumL=99999 # this needs to be set in the C code as well
xPlus=-2./3.
xMinus=4./3.

# ==========================================================================
# C struct setup
# ==========================================================================

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
              ('xPlus', DOUBLE),
              ('xMinus', DOUBLE),
             ]
# define the array of averages that is to be passed to the C routine
paramType=Parameters

# make a class which is a clone of the C struct
# this struct will be passed back with modified arrays
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
# Defining potentials
def Pot(x):
    return 1.+ x*x*(-3.375+x * (1.6875 +x * (2.84766 +(-2.84766+0.711914 * x) * x)))

def Force(x):
    return x * (6.75 + x * (-5.0625 + x * (-11.3906 + (14.2383 - 4.27148 * x) * x)))

def ForcePrime(x):
    return 6.75 + x * (-10.125 + x * (-34.1719 + (56.9531 - 21.3574 * x) * x))

def H(x0,x1,p0):
    dxdt=(x1-x0)*invdt
    return 0.25*(dxdt-Force((x0)))*(dxdt-Force((x0))) + 0.25*(dxdt+Force((x1)))*(dxdt+Force((x1))) + 0.5*p0*p0

def G(x0,x1):
    Fx0=Force((x0))
    Fx1=Force((x1))
    return 0.25*Fx0*Fx0 + 0.25*Fx1*Fx1 + 0.5*(Fx1-Fx0)*(x1-x0)*invdt

def g(xm1,x0,x1):
    g1st=ForcePrime(x0)*Force((x0))
    g2nd=-0.5*invdt*(Force((x1))-2.0*Force((x0))+Force((xm1)))
    g3rd=-0.5*invdt*ForcePrime(x0)*(x1-2.0*x0+xm1)
    return g1st+g2nd+g3rd

# ===============================================
# Gaussian elimination for computing L.x=b where
#   L is tridiagonal matrix
#     mainDiag: 1+2r
#     upper(lower)Diag: -r
#   b is known vector: bVec is input
#   x is unknown vector: returned at the end
def GaussElim(r,bVec):

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

# ===============================================
# Unit Tests
def unitTest(inpath,outpath):
    # print the quadratic variation for every main loop
    quadVar(outpath)
    # verify that the boundary conditions have not changed
    if inpath[0] != outpath[0] and inpath[-1] != outpath[-1]:
        sys.exit("Boundary Conditions Changed! Aborting...")
    
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
params.NumU=NumU
params.NumL=NumL

pathOld=pathType()
pathCur=pathType()
pathNew=pathType()


#matxlist=np.array([0.0 for i in np.arange(1,len(inPath)-1,1)])
#glist=np.array([0.0 for i in np.arange(1,len(inPath)-1,1)])
#noiselist=np.array([0.0 for i in np.arange(1,len(inPath)-1,1)])
#rhs=np.array([0.0 for i in np.arange(1,len(inPath)-1,1)])
outPath=np.array([0.0 for i in np.arange(0,len(inPath),1)])


for loops in range(5):

    for i in np.arange(0,len(inPath),1):
        pathOld[i].pos=inPath[i]
#    # noise vector
    for i in np.arange(1,len(inPath)-1,1):
        pathOld[i].randlist=np.random.normal(0,1)

    c_calcStateSPDE(pathOld,params)

    c_GaussElim(pathOld,pathCur,params)

    for i in np.arange(0,len(inPath),1):
        outPath[i]=pathCur[i].pos
    
    # perform unit tests
#    unitTest(inPath,outPath)

    # save outpath to inpath and loop
    inPath=outPath

if (np.loadtxt("outPathfive.dat") == outPath).all():
    print "SUCCESS!"
else:
    print "FAIL"

np.savetxt("outPathC.dat",outPath)
