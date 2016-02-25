#!/usr/bin/env python3

import numpy as np
import math
import sys

## Sympy generation of the potential and force
## This is hardcoded for this version of the code
#import sympy
#y = sympy.Symbol('y')
#U = sympy.Function('U')
#U = ((8 - 5*y)**8 * (2 + 5*y)**2)/(2**26)
#F = sympy.simplify(-U.diff(y))
#pot = sympy.lambdify(y,U)
#force = sympy.lambdify(y,F)
pot = lambda y : (-5*y + 8)**8*(5*y + 2)**2/67108864
force = lambda y : -125*y*(5*y - 8)**7*(5*y + 2)/33554432

# constant parameters
eps=0.25
dt=0.005
NumB=30001

# basin minima positions
LeftBasin=-2/5.
RightBasin=8/5.

# genStep: function to generate the xnew position
genStep = lambda xold, noise : ( xold + dt*force(xold) + noise )

#===================================
# createTrajectory: main function that performs the simulation.
def createTrajectory():
    
    traj=[0.0 for i in range(NumB)]

    # starting position
    xstart=-2/5.

    # set the starting position 
    xold = xstart

    # prefactor in front of the noise to make v0h
    pref=math.sqrt(2.0*eps*dt)
  
    # accumulators
    acc=0
    rej=0

    #============== MAIN LOOP =====================
    for i in range(NumB):
        # save position to trajectory list
        traj[i]=xold
        
        # generate thermalized noise 
        # noise = v_0 * h = sqrt(2*dt) * sqrt(eps) * xi
        noise = pref*np.random.normal(0.0,1.0)

        # generate the new step (using specific method)
        xnew=genStep(xold, noise)
        
        xold=xnew
        
    return traj


def create_path(broad_frac,broad_frac_precision,ending_precision):
    saved_traj = createTrajectory()
    while True:
        if ( RightBasin-ending_precision < saved_traj[-2] < RightBasin+ending_precision and 
               broad_frac*(1-broad_frac_precision) < sum((np.sign(saved_traj)+1)/2)/len(saved_traj) < broad_frac*(1+broad_frac_precision) ):
            break
        else:
            saved_traj = createTrajectory()
    return saved_traj

np.savetxt("testing.dat",create_path(float(sys.argv[1]),float(sys.argv[2]),math.sqrt(2*eps*dt)))
