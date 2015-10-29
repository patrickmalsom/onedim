# Guided HMC

As a starting point, we will define the potentials in terms of the position $x$

$$ U = U(x,t)$$

$$ U_0 = \frac12 (x(t)-m(t)) A(t) (x(t)-m(t)) - \dot{m} (x(t)-m(t)) $$

And the corresponding stochastic equations are

$$ \frac{dx}{dt} = F + \sqrt{2 \varepsilon} \frac{dW}{dt} $$

$$ \frac{dx}{dt} = -A(t) (x-m(t)) + \frac{dm}{dt} + \sqrt{2 \varepsilon} \frac{dW}{dt} $$

Now performing the hybrid (smart) monte carlo, we have the forward propogation in time

$$ z_1 = z_0 + h v_0 + \frac{h^2}{2} F(\bar{x}) $$
$$ z_1 = z_0 \sqrt{2 \varepsilon ~ \Delta t} - \frac{h^2}{2} A(t+\Delta t /2) \left( \frac{z_1 + z_0}{2} \right) $$

This is an implicit equation for z1, and can be solved iteravely (see SDE methods page). It can also be solved by diagonalizing A (a matrix) and solving for z1 algebraically.

$$ \Delta E = U(z_1+m_1,t_1) - U(z_0+m_0,t_0) + (z_1-z_0) A(t + \Delta t/2) \bar{z} $$

The kinetic term can be understood as the work done (F.d), as the particle feels a constant force between time 0 and 1. The force is given in the equation for z1. now writing the Metropolis adjusted Langevin (MALA) is finalized with the metropolis step itself

$$ \exp \left( -\frac{\Delta E}{\varepsilon} \right) > \eta $$

where $\eta$ is a random number between 0 and 1.

### Solving the implicit equation in one dimension

There are many ways to solve the implicit equation for x1. In one dimension there is a choice that needs to be made on the programming side. One may decide to use the average values of the above formalism (A and m evaluated at the half time step and xhalf=(x0+x1)/2) and then solve for x1, which yields

```
Solve[x1 == x0 + h v0 - h^2/2 (Ahalf ((x1 + x0)/2 - mhalf)), x1]
```

$$ x1 = \frac{ 2 \Delta t ~ A_{1/2} ~ m_{1/2} + 2 \sqrt{2 \Delta t} v0 + 2 x0 - \Delta t ~ A_{1/2} ~ x0}{2 + A_{1/2} \Delta t }$$

alternatively you can use the average A and m evaluated at the end points

```
Solve[x1==x0+hv0-h^2/2 (A1*(x1-m1)+A0*(x0-m0)),x1] $$
```

$$ x1= \frac{ \Delta t ~ A_1 ~ m_1 + \sqrt{2 \Delta t} v0 + \Delta t ~ A_0 (m_0 - x_0) + x_0}{1 + \Delta t ~ A_1}$$

When dealing with multiple dimensions, the implicit methods above will need to be solved either by some iterative method or by diagonalizing the A matrix for the step. This will obviously need some more development. As a side note, when doing this in multiple dimensions, you can use leap frog (velocity at the half step) which allows the equation to be solved directly (but adds to the error).


# Example 1D potential

At this point it is ilustrative to show an example one dimensional potental with the corresponding $m(t)$ and $A(t)$.

### insert image of potential 

The potential $U(x)= (x^2 - 1 ) ^2 $

### insert m(t) image

The mean m(t) is defined to be $m^-$ in teh starting well (t=0), and $m^+$ in the final well (t=T). This function is constant for almost half of the path at which point it quickly transitions to the other well and resumes its constant behaviour.

### insert A(t) image

The quadratic term in the potential examansion is a constant balue $A_0$, while sittin in either well (for this specific potential) and changes rapidly in teh transition region.

The values $m^+$, $m^-$ and $A_0$ can hopefully be found using the same analysis performed in the warwick notebook. It is yet to be seen what values the forward integration will yield (free energy, path space, other...). This m)t) and A(t) will "drive" the trasition and translate it into a stochastic control problem. Look at the literature on driven or steered molecular/stochastic dynamics.


 * we want an "optimized" A(t) and m(t). What does that mean?
 * generate a "path" using equation 1.
 * this requires an initial m an A and then need to iterate to find improved m and A.
 * the sequence of configurations (multiple paths) are not necessarily driven by the dynamics.
 * we want the distribution with m and A to appr oximate the boltzmann distribution witht the true U(x).
 * optimization means minimizing the KL distance between the two distributions.

$$Z_B = \int \mathcal{D} z \int dt \exp \left( - \frac{U(z+m)}{ \varepsilon} \right) $$

 * $\mathcal{D} z$ corresponds tot the sequences or "paths".
 * $dt$ corresponds tot eh time along the sequence.
 * only in the continuum limit does this correspond to the true dynamics.

$$Z_0 = \int \mathcal{D} z \int dt \exp \left( - \frac{U_0(x,t)}{ \varepsilon} \right) $$

# Jensens Inequality

$$ \frac{ \mathcal{Z}_B }{ \mathcal{Z}_0} = \mathbb{E}_0 \exp \left( - \frac{ (U-U_0)}{\varepsilon} \right) $$

$$ \frac{ \mathcal{Z}_B }{ \mathcal{Z}_0} \geq  \exp \left( - \frac{1}{\varepsilon} \mathbb{E}_0 (U-U_0) \right) $$

$$ \ln \left( \frac{ \mathcal{Z}_B }{ \mathcal{Z}_0} \right) = \ln (\mathcal{Z}_B) - \ln (\mathcal{Z}_0 ) $$

$$ - \varepsilon \ln (\mathcal{Z}_B ) \leq -\varepsilon \ln (\mathcal{Z}_0) + \mathbb{E}_0(U-U_0) $$

$$ F_B \leq F_0 + \mathbb{E}_0 (U-U_0) $$

Note that $F_0$ is independent of m and $F_B$ is independent of A and m.

# KL Distance

$$ \varepsilon D_{KL} = F_0 - F_B + \mathbb{E} (U-U_0) \geq 0 $$

The KL distance is very similar to Jensens Inequality

$$D_{KL} = \ln \mathcal{Z}_B - \ln \mathcal{Z}_0 + \frac{1}{\varepsilon} \mathbb{E}_0 (U-U_0) $$

using the substitution $\Delta U = U - U_0 $

$$ D_{KL} = \mathbb{E}_0 \left( \frac{\Delta U}{\varepsilon} \right) + \ln \left( \mathbb{E}_0 \left( \exp \left( - \frac{\Delta U}{\varepsilon} \right) \right) \right) $$

# Optimization of A

To review, we have 
$$U = U(x)$$ 
$$U_0 = \frac12 (x-m)A(x-m) - \dot{m}x$$

and the Kullback-Leibler distance defined to be

$$ D_{KL} = \frac{1}{\varepsilon} \mathbb{E}_0 (U-U_0) + \ln \left( \frac{Z}{Z_0} \right) $$

$$ \frac{ \partial D_{KL} }{\partial \alpha} = \frac{1}{\varepsilon} \mathbb{E}_0(U-U_0) + Z_0 \frac{\partial}{\partial \alpha} \frac{1}{Z_0} $$

$$ \frac{\partial}{\partial \alpha} \frac{1}{Z_0} $$

# Problems with the guided routine

There seems to be many problems with the guided method. All of these problems stem from the causal nature of the dynamics. In order to explain this I will talk about two sets of objects, controling objects and response objects. The controlling objects in this algorithm are m(t) and A(t), acting as a force to push the particle from one well to another. These objects can be set to be any function and a set of responses are calculated from them.

A response is any object calculated from the control, and all of these responses have causality built into them, which can be thought of as a lag in time from the control. Examples of the response are xbar and all of the KL derivatives. 

There is a fine line to walk when attempting to force transition to occur by imposing some mathematical framework. At the end of this work, we have come to realize that the minimum of the KL distance for a particle sitting in a well is for the particle to sit in the well forever. This might seem obvious at first glance. This behaviour is shown to occur for the above algorithm, when the control functions being forced to make the transition at longer and longer times. 

There were many ideas about how to deal with the disparate nature of the control and response, none of which have worked so far. One idea was to shift the response backward in time to match up with the control functions, which failed to converge when using gradient descent to minimize the KL distance. 
