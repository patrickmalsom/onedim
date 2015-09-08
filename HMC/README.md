# HMC algorithm for a 1D potential

This code runs the HMC algorithm for a particle in a 1D potential.
The algorithm is implimented for the following methods:

|  method  | description |
| -------- | ----------- |
| ito      | Ito/Girsanov form of HMC (continuous time limit)                |
| midpt    | Finite time HMC using the mid point integrator                  |
| leapfrog | Finite time HMC using the LeapFrog (velocity Verlet) integrator |

## Overview of the code

This code uses python's ctypes to exploit the speed of a compiled C library.

## Compiling the C library

Before running this code, the C library must be compiled using the `make` command.

## Running the simulation

`onedimHMC.py` is the main executable used to run the simulation.
to print the help run the following

```
./onedimHMC.py -h
```

### Adding a new potential

Including another potential to the library is very simple

1. generate the forces for the potential
    - Force
    - ForcePrime (dF/dx)
    - ForceDoublePrime (d^2F/dx^2)
2. add these calls to a c file in the `potential_defns` directory
3. include the following header files as includes

``` c
#include "../onedimHMC_struct.h"
#include "../onedimHMC_forces.h"
#include "../onedimHMC.c"
```

4. run the program with the `--potential` flag corresponding to the name of the file in the `potential_defns` directory

```
./onedimHMC.py --method ito --potential fatter_skinny -i input_paths/fatterSkinny-T0p25-dt0p005-Nb30k-healed.dat -T 0.25 --deltat 0.005 --deltatau 0.0000001 --Num 30001 --HMC 10 --MD 1000
```
