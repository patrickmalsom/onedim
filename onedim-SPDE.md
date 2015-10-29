# SPDE

The differential equation is given as

$$\frac{\partial x}{\partial t}=\frac{\partial\^2 x}{\partial u\^2} - \frac{\partial G}{\partial x}+ \sqrt{2 \varepsilon} \frac{dW }{dt}$$

Where $G$ is the path potential, defined to be

$$ G=\frac{1}{2} \left( |F|\^2 -2 \varepsilon \nabla\^2 V \right) $$

Crank nicholson discretization yields

$$ \frac{z_i-x_i}{\Delta t} = \frac12 \left( \frac{x_{i-1} - 2 x_i + x_{i+1}}{\Delta u\^2} + \frac{z_{i-1}-2 z_i+z_{i+1}}{\Delta u\^2} \right)-\nabla \left(G x_i \right) + \xi_i \sqrt{\frac{4 \varepsilon}{\Delta t \Delta u\^2}} $$

In the above equation z is future time and x is current time, $\xi$ is a gaussian random number and the incriment i is the path time (du). Now if we separate future time from current time (z from x) and substitute $r=\frac{\Delta t}{2 \Delta u\^2}$ we find

$$z_i - r \left( z_{i-1}+z_{i+1}-2 z_i \right) = \Delta t \left( - \left( G x_i \right) \right) + r \left( x_{i-1}+x_{i+1} \right) + (1-2 r) x_i+ \xi_i \sqrt{4 \varepsilon \frac{\Delta t}{\Delta u\^2}} $$

This is a matrix equation of the form $(I+M)z=b$ where M is tridiagonal and looks like

$$
\left( \begin{array}{ccccc}
 2 & -1 & 0 & 0 & 0 \\\\
 -1 & 2 & -1 & 0 & 0 \\\\
 0 & -1 & 2 & -1 & 0 \\\\
 0 & 0 & -1 & 2 & \ddots \\\\
 0 & 0 & 0 & \ddots & \ddots
\end{array} \right) 
$$

The routine to run the SPDE goes as follows 

 1. Generate a generic starting path (hopefully something close to the physical path for reasonable convergence times)
 2. use the above descrete equation to solve the matrix equation and thus solve for the future positions (z). This involves gaussian elimination on the lower off diagonal and then back substitution to find the final postions. You must be careful with the boundary conditions of the matrix. We have ommited some information by making the matrix square.
 3. Save the newly found z configuration to the x postions and repeat.

## generation of a path

The next step is to generate the starting path for the SPDE. This path can be anything but the endpoint will never move due to boundary conditions. This means that the path should probably start and end in local or global minima. 

```
(* starting guess for x *)

leftxMin = -1;(*position of the left minimum in the energy*)

rightxMin = 1;(*position of the right minimum in the energy*)

steepx = 10 du;
posx = Table[(rightxMin - 
       leftxMin) (Tanh[steepx (x/du - (numBead - 1)/2)]/2 Tanh[
         steepx*(numBead - 1)/2] + 0.5) + leftxMin, {x, 
    0, (numBead - 1)*du, du}];

ListLinePlot[posx, Frame -> True, Axes -> False]
```


```
(* starting guess for y *)

leftyMin = 0;(*position of the left minimum in the energy*)

rightyMin = 0;(*position of the right minimum in the energy*)

steepy = 10 du;
posy = Table[(rightyMin - 
       leftyMin) (Tanh[steepy (y/du - (numBead - 1)/2)]/2 Tanh[
         steepy*(numBead - 1)/2] + 0.5) + leftyMin, {y, 
    0, (numBead - 1)*du, du}];

ListLinePlot[posy, Frame -> True, Axes -> False]
```

```
SetAttributes[SPDE, HoldAll];
SPDE[posListx_, posListy_, du_, Temp_, dt_, numBead_] :=
 
 Module[{pref = Sqrt[4 Temp dt/du], r = dt/(2 du^2), A, C2},
  (* X positions *)
  
  A = Table[{ 1 + 2 r, 
     posListx[[j]] + 
      r   (posListx[[j + 1]] - 2 posListx[[j]] + posListx[[j - 1]]) - 
      dt*xdelG[posListx[[j]], posListy[[j]], Temp] + 
      pref RandomReal[NormalDistribution[]]}, {j, 2, numBead - 1}];
  A[[1, 2]] += r posListx[[1]];
  A[[Length[A], 2]] += r posListx[[ Length[posListx]]];
  Do[ A[[i, 1]] -= r^2/A[[i - 1, 1]];
   A[[i, 2]] += r A[[i - 1, 2]]/A[[i - 1, 1]];
   , {i, 2, numBead - 2}];
  C2 = Table[0, {i, 1, numBead - 2, 1}];
  C2[[numBead - 2]] = A[[numBead - 2, 2]]/A[[numBead - 2, 1]];
  Do[C2[[i - 1]] = (A[[i - 1, 2]] + r C2[[i]])/A[[i - 1, 1]];
   , {i, numBead - 2, 2, -1}];
  PrependTo[C2, First[posListx]];
  AppendTo[C2, Last[posListx]];
  posListx = C2;
  
  (* Y positions *)
  
  A = Table[{ 1 + 2 r, 
     posListy[[j]] + 
      r   (posListy[[j + 1]] - 2 posListy[[j]] + posListy[[j - 1]]) - 
      dt*ydelG[posListx[[j]], posListy[[j]], Temp] + 
      pref RandomReal[NormalDistribution[]]}, {j, 2, numBead - 1}];
  A[[1, 2]] += r posListy[[1]];
  A[[Length[A], 2]] += r posListy[[ Length[posListy]]];
  Do[ A[[i, 1]] -= r^2/A[[i - 1, 1]];
   A[[i, 2]] += r A[[i - 1, 2]]/A[[i - 1, 1]];
   , {i, 2, numBead - 2}];
  C2 = Table[0, {i, 1, numBead - 2, 1}];
  C2[[numBead - 2]] = A[[numBead - 2, 2]]/A[[numBead - 2, 1]];
  Do[C2[[i - 1]] = (A[[i - 1, 2]] + r C2[[i]])/A[[i - 1, 1]];
   , {i, numBead - 2, 2, -1}];
  PrependTo[C2, First[posListy]];
  AppendTo[C2, Last[posListy]];
  posListy = C2;
  ]
```

```
(* Perform the SPDE simulation *)
AbsoluteTiming[
 Do[SPDE[posx, posy, du, Temp, dt, numBead], {i, 1, 100}]]
```
