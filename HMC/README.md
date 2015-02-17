# HMC algorithm for a 1D potential

This code runs the HMC algorithm for a particle in a 1D potential.

It includes the following methods

|  method  | description |
| ======== |:=========== |
| ito      | Ito/Girsanov form of HMC (continuous time limit) |
| midpt    | Finite time HMC using the mid point integrator |
| leapfrog | Finite time HMC using the LeapFrog (velocity Verlet) integrator |

### Adding new potentials

Including another potential to the library is very simple

1. Make a new `potential.c` file which declares the forces
2. add the new potential to the `Makefile`
3. add the new potential name to argparse in `onedimHMC.py`
4. add a new try statement to import the c library
