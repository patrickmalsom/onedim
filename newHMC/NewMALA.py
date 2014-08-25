#!/usr/bin/env python2.7

# MALA routine

# ===============================================
# Import libraries
import numpy as np
import math
import sys


# ===============================================
# constants/parameters
deltat=0.005
invdt=1/deltat
eps=0.15
deltatau=0.00001
noisePref=math.sqrt(4.0*eps*deltatau*invdt)
np.random.seed(101)

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
    return 0.25*(dxdt-Force(x0))*(dxdt-Force(x0)) + 0.25*(dxdt+Force(x1))*(dxdt+Force(x1)) + 0.5*p0*p0

def G(x0,x1):
    Fx0=Force(x0)
    Fx1=Force(x1)
    return 0.25*Fx0*Fx0 + 0.25*Fx1*Fx1 + 0.5*(Fx1-Fx0)*(x1-x0)*invdt

def g(xm1,x0,x1):
    g1st=ForcePrime(x0)*Force(x0)
    g2nd=-0.5*invdt*(Force(x1)-2.0*Force(x0)+Force(xm1))
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

    tempVec=np.array(bVec)
    Num=len(tempVec)
    
    al=[-r for i in range(Num-1)]
    am=[(1+2*r) for i in range(Num)]
    au=[-r for i in range(Num-1)]


    # Gaussian Elimination
    for i in range(Num-1):
        temp=-al[i]/am[i]
        #al[i]+=am[i]*temp
        am[i+1]+=au[i]*temp
        tempVec[i+1]+=tempVec[i]*temp

    # Back substitution
    for i in range(Num-1):
        temp=-au[Num-i-2]/am[Num-i-1]
        #au[Num-i-2]+=am[Num-i-1]*temp
        tempVec[Num-i-2]+=tempVec[Num-i-1]*temp
    
    # Divide by main diagonal
    for i in range(Num):
        tempVec[i]=tempVec[i]/am[i]
        #am[i]=am[i]/am[i]
    
    # Print the result
    return tempVec

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

inPath=np.loadtxt("SDEPath-fat-skinny.dat")

for loops in range(100):

    # x vector
    matxlist=[0.0 for i in np.arange(1,len(inPath)-1,1)]
    r=deltatau*0.5*invdt*invdt
    for i in np.arange(1,len(inPath)-1,1):
        matxlist[i-1]=r*inPath[i-1] + (1.-2.*r)*inPath[i] + r*inPath[i+1]

    # g vector
    glist=[0.0 for i in np.arange(1,len(inPath)-1,1)]
    for i in np.arange(1,len(inPath)-1,1):
        glist[i-1]= deltatau * g(inPath[i-1],inPath[i],inPath[i+1])
    
    # noise vector
    noiselist=[0.0 for i in np.arange(1,len(inPath)-1,1)]
    for i in np.arange(1,len(inPath)-1,1):
        noiselist[i-1]= noisePref*np.random.normal(0,1)
    
    # make full rhs vector
    rhs=[0.0 for i in np.arange(1,len(inPath)-1,1)]
    for i in np.arange(0,len(inPath)-2,1):
        rhs[i]= matxlist[i]-glist[i]+noiselist[i]
    # adding the BC from the LHS matrix equation
    rhs[0]+=r*inPath[0]
    rhs[-1]+=r*inPath[-1]

    sol=GaussElim(r,rhs)

    # add BC's back into solution
    outPath=np.append(np.insert(sol,0,inPath[0]),inPath[-1])
    
    # perform unit tests
    unitTest(inPath,outPath)

    # save outpath to inpath and loop
    inPath=outPath
