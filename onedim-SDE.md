# Forward Integration of the SDE in 1D

This is the wiki page for a set of programs designed to simulate a particle in a one dimensional potential. The program uses 3 separate schemes for the integration step

 * Euler (Ito)
 * Mid Point
 * Rotation

furthermore, the code can also perform a Metropolis-Hastings Monte-Carlo step if required. All of these options will be outlined below

### Metropolis Hastings
A guarantee that the underlying probability distribution is Gaussian can be obtained by using a Metropolis-Hastings test to force the probability to be $\mathcal{P} \propto exp(-U/\varepsilon)$. The Metropolis condition in pseudocode is something like:

``` python
if (exp(-DeltaE/(2*eps))>rand):
  x0=x1
  acc+=1
else
  x0=x0
  rej+=1
```

The Metropolis adjusted sampling requires the calculation of the change in energy when performing the move from $x_0$ to $x_1$ ($\Delta E$). This energy change (error) is written simply as

$$ \Delta E = U(x_1) - U(x_0) + \frac{1}{2} v_1\^2 - \frac{1}{2} v_0\^2 $$

$x_0$ is known (or guessed for the first step) and $v_0$ is chosen from a random Gaussian distribution, and $x_1$ is calculated using the current position and velocity. The integration methods of interest to this discussion are not aware of the future velocities ($v_1$) at all, as the future velocity will be chosen as the current velocity when performing the next step. For this reason, the future velocity must be eliminated from the calculations for the energy error for the current step. The energy error associated with a particular move is dependent on the employed integration method, as is shown seperately for each integration method in the sections below.

# Integration Methods

----------------------------------------------
### Euler (Ito)

In the over-damped limit, the equation of motion of the particle is given by 

$$ \frac{dx}{dt} = F(x) + Noise $$

The generation of the first step is performed by picking a starting $x_0$ and evolving it forward in time. 

$$ x_1 = x_0 + \Delta t F(x_0) + \sqrt{2 \varepsilon ~ \Delta t } \xi_0 $$

where subscript 0 denotes the current step and subscript 1 denotes the future (generated) step, and $\xi$ is a Gaussian distributed random number. Note that the future position is defined entirely by the current position and velocity.

Scaling the time and the velocity to be $h=\sqrt{2 ~ \Delta t}$ and $v_0=\sqrt{\varepsilon} \xi_0$ respeciively, we find the familiar form for the equation of motion of the particle

$$ x_1 = x_0 + h ~ v_0 + \frac12 h\^2 F(x_0) $$

Performing the leap frog evaluation of the velocity (evaluate at the midpoint) yields $v_1$

$$v_1 = v_0 + \frac{h}{2} F \left( \frac{x_1 + x_0}{2} \right) $$

In order to perform the MH step we need the total change in energy when moving from $x_0$ to $x_1$.

$$\Delta E = U(x_1) - U(x_0) + \Delta KE$$

The kinetic energy term can be written in terms of the current and future postions as is shown below

$$ \begin{aligned} 
\Delta KE &= \frac12 v_1\^2 - \frac12 v_0\^2 \newline
&= \frac12 \left( v_0 + \frac{h}{2} (F_1 + F_0) \right)\^2 - \frac12 V_0\^2 \newline
&= \frac{h}{2} v_0  ( F_1+F_0 ) + \frac{h\^2}{8} (F_1+F_0)\^2 \newline
&= \frac12 (x_1-x_0)(F_1 + F_0) - \frac14  h\^2 F_0 (F_1+F_0) + \frac{h\^2}{8} (F_1+F_0)\^2 \newline
&= \frac12 (x_1-x_0)(F_1 + F_0) + \frac{h\^2}{8} (F_1+F_0)(-2F_0 + F_1 + F_0) \newline
&= \frac12 (x_1-x_0)(F_1 + F_0) + \frac{h\^2}{8} (F_1\^2+F_0\^2) \end{aligned} 
$$

----------------------------------------------
### Mid Point

The Mid-Point method is very similar to the to the Euler-Ito iteration outlined above, other than the point at which the forces are evaluated. Starting with the discrete equation

$$ x_1 = x_0 + v_0 h + \frac12 h\^2 F \left( \frac{x_1+x_2}{2} \right) $$
$$ v_0 = \sqrt{\varepsilon} \xi_0 $$
$$v_1 = v_0 + h F \left( \frac{x_1+x_2}{2} \right) $$

These give the following velocities

$$v_0 h = x_1 - x_0 - \frac12 h\^2 F \left( \frac{x_1+x_2}{2} \right) $$
$$v_1 h = x_1 - x_0 + \frac12 h\^2 F \left( \frac{x_1+x_2}{2} \right) $$

The sum and the difference between the velocities are useful when calculating the kinetic energy change

$$ v_1 + v_0 = \frac{2}{h} (x_1 - x_0) $$
$$ v_1 - v_0 = h F \left( \frac{x_1+x_2}{2} \right) $$

now the kinetic energy is given as
$$ \begin{aligned} 
\Delta KE &= \frac12 v_1\^2 - \frac12 v_0\^2 \newline
&= \frac12 (v_1 + v_0) (v_1 - v_0) \newline
&=  (x_1 - x_0) F \left( \frac{x_1 + x_0}{2} \right) \end{aligned} 
$$

#### Notes on Mid-Point Method

This method requires the evaluation of the force at the midpoint before determining the final position $x_1$. Normally, the force is complicated enough that separation of $x_1$ is very difficult or impossible. To remedy this, I have use an iterative scheme to generate the new $x_1$. The first step is to guess $x_1$ using the iterative scheme and then use this new guess to generate the next iteration of the future position. Although 3 or 4 of these iterations are usually sufficient to find $x_1$ to machine precision, this scheme is not gaurenteed to converge and should be watched closely. Pseudocode for the schme follows

``` python
x1 = x0 + v0*h + 0.5*h**2*F(x0)
for i in range(4):
  guessx = x0 + v0*h + 0.5*h**2*F((x0+x1)/2)
  x1 = guessx
```

----------------------------------------------
### Rotation
Another way to proceed is to simply solve the differential equation for the case of a constant fluctuation. The ODE is

$$ x''(t) = F - (x(t) - \bar{x} ) A $$

The homogeneous solution

$$ x'' + A x = 0 $$
$$ x_h = c_1 sin(\sqrt{A} t )  + c_2 cos(\sqrt{A} t) $$

The particular solution ($\beta = F + \bar{x} A$)

$$ x'' + A x = \beta $$
now guess that the solution is of the form
$$ x_p = \gamma \rightarrow A \gamma = \beta \rightarrow \gamma = \frac{\beta}{A} $$
$$ x_p = \frac{F}{A} + \bar{x} $$

Full solution

$$ x = x_h + x_p $$
$$ x = c_1 sin(\sqrt{A} t ) + c_2 cos( \sqrt{A} t) + \frac{F}{A} + \bar{x} $$

Imposing the boundary conditions $x(0) = x0$ and $\dot{x}(0) = v_0$

$$ x(0) = c_2 + \frac{F}{A} + \bar{x} = x_0 $$
$$ \dot{x}(0) = \sqrt{A} c_1 = v_0 $$

Finally arriving at the full solution

$$ x(t) = x_0 cos(\sqrt{A} t ) + \frac{v_0}{\sqrt{A}} sin( \sqrt{A} t ) + \left( \frac{F}{A} + \bar{x} \right) \left( 1 - cos( \sqrt{A} t) \right) $$

The discretization can use either the midpoint or the euler method. Evaluation at the midpoint 

$$ x(h) = x_1 = \left( \frac{F(\bar{x})}{A(\bar{x})} + \bar{x} \right) \left( 1 - cos \left(\sqrt{2 A(\bar{x}) \Delta t} \right) \right) + \frac{\sqrt{\varepsilon} \xi_0 }{\sqrt{A(\bar{x})}} sin \left( \sqrt{2 A(\bar{x}) \Delta t} \right) + x_0 cos \left( \sqrt{2 A(\bar{x}) \Delta t} \right) $$

where $\bar{x} = (x_1+x_0)/2 $. This will have to be iterated on similarly to the midpoint case.


# Observations and Thoughts

* The quadratic variation  is **very** important to look at. For the non-Metropolis adjusted method, the quadratic variation goes as $V''(x)$. This implies that the Metropolis adjusted case will show many rejections in the region that the value of $|V''(x)|$ is large and $V(x)$ is small (compared to their minimum values).

* The Metropolis step is useless at small $\Delta t$ (as expected) as the acceptance rate is close to unit, and thus Metropolis Hastings is never even relied upon in this regime. When $\Delta t$ is large, the non-Metropolis method **tends to diverge** as the particle will be forced up the side of a well (and cannot be rejected) where the subsequent move will have a large force and the resulting move will overshoot the well entirely.

* For any of the midpoint methods, it is important to understand the computational cost of iterating on a complex iterator. The error in energy is smaller for the midpoint and rotation methods for most cases, but the computational cost of computing $x_1$ iteratively counterbalances this savings. Furthermore, the rotation method is **exact** for the OU (quadratic well) case but the cost to calculate $x_1$ is very expensive, as there are sqrt and trig functions in the forward integration.
