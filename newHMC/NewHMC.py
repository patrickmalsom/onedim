#!/usr/bin/env python2.7

import numpy as np
import math

# =================================================================
# Constants
deltat=0.001
invdt=1/deltat
eps=0.15
deltatau=0.00001
noisePref=math.sqrt(4.0*eps*deltatau*invdt)
r=deltatau*0.5*invdt*invdt

# =================================================================
# ## Read in Path
inPath=np.loadtxt("SDEPath-fat-skinny.dat")
Num=len(inPath)

# =================================================================
# seed the RNG
np.random.seed(100)

# =================================================================
# Function declarations
def Pot(x):
    return 1.+ x*x*(-3.375+x * (1.6875 +x * (2.84766 +(-2.84766+0.711914 * x) * x)))

def Force(x):
    return x * (6.75 + x * (-5.0625 + x * (-11.3906 + (14.2383 - 4.27148 * x) * x)))

def ForcePrime(x):
    return 6.75 + x * (-10.125 + x * (-34.1719 + (56.9531 - 21.3574 * x) * x))

# ## Effective Hamiltonian
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

def quadVar(poslist):
    print "quad var: " + str(sum([(poslist[i]-poslist[i+1])**2 for i in range(len(poslist)-1)])/(2.*eps*deltat*(len(poslist)-1)))


def GaussElim(bVec):
# Solves Mx=b where
#  M is tridiagonal matrix
#    mainDiag: Length Num-2 (BC's)
#    upper(lower)Diag: Length Num-3
#  x is unknown and returned at the end
#  b is known vector

    # size of the matrix
    mainLen=len(bVec)

    # make the main/lower/upper diagonal lists
    mainDiag=[1+2*r for i in range(mainLen)]
    lowerDiag=[-r for i in range(mainLen-1)]
    upperDiag=[-r for i in range(mainLen-1)]

    # Gaussian Elimination down the lower diag
    for i in range(mainLen-1):
        temp=-lowerDiag[i]/mainDiag[i]
        #lowerDiag[i]+=mainDiag[i]*temp
        mainDiag[i+1]+=upperDiag[i]*temp
        bVec[i+1]+=bVec[i]*temp

    # Gaussian Elimination up the upper diag
    for i in range(mainLen-1):
        temp=-upperDiag[mainLen-i-2]/mainDiag[mainLen-i-1]
        #upperDiag[mainLen-i-2]+=mainDiag[mainLen-i-1]*temp
        bVec[mainLen-i-2]+=bVec[mainLen-i-1]*temp

    # Divide by main diagonal
    for i in range(mainLen):
        bVec[i]=bVec[i]/mainDiag[i]
        #mainDiag[i]=mainDiag[i]/mainDiag[i]

    # return the result
    return bVec

def genRHS(inputpath):
    # (I+hL)x0
    matxlist=0.0
    # g_n
    glist=0.0
    # noise
    noiselist=0.0
    knownVec=[0.0 for i in np.arange(1,len(inPath)-1,1)]

    # calculate the nth element of the known vector on the RHS
    for i in np.arange(1,len(inPath)-1,1):
        matxlist=r*inPath[i-1] + (1-2.*r)*inPath[i] + r*inPath[i+1]
        glist= deltatau * g(inPath[i-1],inPath[i],inPath[i+1])
        noiselist= noisePref*np.random.normal(1,1)
        knownVec[i-1]= matxlist-glist+noiselist

    # adding the BC from the LHS matrix equation
    knownVec[0]+=r*inPath[0]
    knownVec[-1]+=r*inPath[-1]

    return knownVec

# =================================================================
# generate the known RHS vector
rhs=genRHS(inPath)

# perform the gaussian elimination to find the new path x^1
xnewPath=GaussElim(rhs)

# print the quadratic variation
quadVar(xnewPath)
# save the output path to file
np.savetxt("NewHMC-out.dat",xnewPath)


# =================================================================
# ## Molecular Dynamics
# 
# Rather than solve for the next steps explicitly in terms of the previous two steps, as we did for the path space implementation, we are going to solve for the velocities and then use the same algorithm as the MALA above.
# This involves solving explicitly for the new velocities from the old noise list ($h~p_0$).
# 
# $$p^{(1)} = p^{(0)} + \sqrt{ \frac{\Delta \tau}{2}} \left( L x^{(1)} + L x^{(0)} \right) + \sqrt{ \frac{\Delta \tau}{2}} \left( g^{(1)}-g^{(0)} \right) $$
# 
# After finding the velocities for the next step we can return to the top to perform the same *MALA* without choosing random velocities, thus producing the MD step.
# 
# ## Calculation of $L x$
# 
# The matrix L will need to be multiplied by a vector many times in this algorithm.
# I will therefore make a function to do so.
# 
# ## Fourier Transform
# 
# * probably should use $2^n +1$ points, allows easy use of FFT
# 
# $\Delta x$ transform should show high frequency behaviour
#  
# $$\sqrt{ \frac{\Delta t}{ 2 \epsilon}} \left( \frac{\Delta x}{\Delta t} - F \right) = \frac{A_{\omega}}{\omega}$$
# 
# This should therefore show the $1/f$ noise that is expected of the process. $\sum A_{\omega}^2 = N$ and  $\sum A_{\omega} = 0$
# 
# The question is: will this method heal the high frequency stuff to be correct if the input guess is crap?


