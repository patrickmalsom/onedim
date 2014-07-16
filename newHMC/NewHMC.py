#!/usr/bin/env python2.7

import numpy as np
import math
from pylab import *

# Constants
deltat=0.001
invdt=1/deltat
eps=0.15
deltatau=0.00001
noisePref=math.sqrt(4.0*eps*deltatau*invdt)

# ## Read in Path
inPath=loadtxt("SDEPath-fat-skinny.dat")
Num=len(inPath)

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


# Solves Mx=b where
#  M is tridiagonal matrix
#    mainDiag: Length Num-2 (BC's)
#    upper(lower)Diag: Length Num-3
#  x is unknown and returned at the end
#  b is known vector
def GaussElim(mainDiag,lowerDiag,upperDiag,bVec):
    # size of the matrix
    mainLen=len(mainDiag)

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


# ## generate x0 vector on RHS
matxlist=[0.0 for i in arange(1,len(inPath)-1,1)]
r=deltatau*0.5*invdt*invdt
for i in arange(1,len(inPath)-1,1):
    matxlist[i-1]=r*inPath[i-1] + (1-2.*r)*inPath[i] + r*inPath[i+1]

# ## generate g list
glist=[0.0 for i in arange(1,len(inPath)-1,1)]
for i in arange(1,len(inPath)-1,1):
    glist[i-1]= deltatau * g(inPath[i-1],inPath[i],inPath[i+1])


# ## generate noise list
noiselist=[0.0 for i in arange(1,len(inPath)-1,1)]
np.random.seed(100)
for i in arange(1,len(inPath)-1,1):
    noiselist[i-1]= noisePref*np.random.normal(1,1)


# ## make full RHS vector
rhs=[0.0 for i in arange(1,len(inPath)-1,1)]
for i in arange(0,len(inPath)-2,1):
    rhs[i]= matxlist[i]-glist[i]+noiselist[i]
# adding the BC from the LHS matrix equation
rhs[0]+=r*inPath[0]
rhs[-1]+=r*inPath[-1]

# perform the gaussian elimination to find the new path x^1
md=[1+2*r for i in range(len(inPath)-2)]
ld=[-r for i in range(len(inPath)-3)]
ud=[-r for i in range(len(inPath)-3)]
xnewPath=GaussElim(md,ld,ud,rhs)

# print the quadratic variation
quadVar(xnewPath)

# save the output path to file
savetxt("NewHMC-out.dat",xnewPath)


# ## Molecular Dynamics
# 
# Rather than solve for the next steps explicitly in terms of the previous two steps, as we did for the path space implementation, we are going to solve for the velocities and then use the same algorithm as the MALA above.
# This involves solving explicitly for the new velocities from the old noise list ($h~p_0$).
# 
# $$p^{(1)} = p^{(0)} + \sqrt{ \frac{\Delta \tau}{2}} \left( L x^{(1)} + L x^{(0)} \right) + \sqrt{ \frac{\Delta \tau}{2}} \left( g^{(1)}-g^{(0)} \right) $$
# 
# After finding the velocities for the next step we can return to the top to perform the same *MALA* without choosing random velocities, thus producing the MD step.
# 
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


