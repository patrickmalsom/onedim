
# coding: utf-8

# In[79]:

get_ipython().magic('matplotlib inline')
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
import sympy

from matplotlib import rc_file
rc_file('/home/patrick/.matplotlib/matplotlibrc')


# In[94]:

y = sympy.Symbol('y')
U = sympy.Function('U')
U = ((8 - 5*y)**8 * (2 + 5*y)**2)/(2**26)
F = sympy.simplify(-U.diff(y))


# In[95]:

pot = sympy.lambdify(y,U)
force = sympy.lambdify(y,F)


# In[154]:

eps=0.25
dt=0.005
NumB=30001

# basin minima
LeftBasin=-2/5.
RightBasin=8/5.


# In[137]:

# genStep: function to generate the xnew position
genStep = lambda xold, noise : ( xold + dt*force(xold) + noise )

#===================================
# calcEnergy: function that calculates the energy error associated
#    with the integration forward. 
#def calcEnergy(xold, xnew):
#        # evaluate the force at current and next step for later use
#        Fx1, Fx0 = force(xnew), force(xold)
#        # return the energy error. see notes for explanation of algebra.
#        return 0.5*(xnew-xold)*(Fx1+Fx0)+dt*0.25*(Fx1*Fx1-Fx0*Fx0) + (pot(xnew) - pot(xold))


# In[138]:

genStep(0,.1)


# In[139]:

plt.plot([x for x in np.arange(-.7,3,0.01)],[pot(x) for x in np.arange(-.7,3,0.01)])
plt.ylim(-0.1,2)
plt.title(r"Potential  $U = [(8 - 5y)^8  (2 + 5y)^2]/2^{26}$")
plt.show()

plt.plot([x for x in np.arange(-.7,3,0.01)],[force(x) for x in np.arange(-.7,3,0.01)])
plt.ylim(-5,3)
plt.title("Force")
plt.show()


# In[155]:

#===================================
# createTrajectory: main function that performs the simulation.
#    prints results to stdout when finished
def createTrajectory():
    
#    if traj_len < NumB:
#        print("need longer trajectory length")
#        break
    
    traj=[0.0 for i in range(NumB)]

    # starting position
    xstart=-2/5.

    # set the starting position 
    xold = xstart

    # prefactor in front of the noise to make v0h
    pref=math.sqrt(2.0*eps*dt)
  
    # accumulators
    trans=0
    acc=0
    rej=0
    


    #============== MAIN LOOP =====================
    for i in range(NumB):
        # save position to trajectory list
        traj[i]=xold
        
        # generate thermalized noise to use when generating the new step 
        #   noise = v_0 * h = sqrt(2*dt) * sqrt(eps) * xi
        noise = pref*np.random.normal(0.0,1.0)

        
        # generate the new step (using specific method)
        xnew=genStep(xold, noise)
        
        xold=xnew
        
        
    return traj


# In[218]:

def create_path(broad_frac,broad_frac_precision,ending_precision):
    saved_traj = createTrajectory()
    while True:
        if ( RightBasin-ending_precision < saved_traj[-2] < RightBasin+ending_precision and 
               broad_frac*(1-broad_frac_precision) < sum((np.sign(saved_traj)+1)/2)/len(saved_traj) < broad_frac*(1+broad_frac_precision) ):
            break
        else:
            saved_traj = createTrajectory()
    return saved_traj


# In[219]:

sum((np.sign(saved_traj)+1)/2)/len(saved_traj)


# In[224]:

saved_traj=create_path(0.7,0.05,math.sqrt(2*eps*dt))
plt.plot(saved_traj)
plt.show()


# In[228]:

saved_traj[-2]


# In[ ]:



