# HMC algorithm for a 1D potential

This code runs the HMC algorithm for a particle in a 1D potential.

It includes the following methods

|  method  | description |
| -------- | ----------- |
| ito      | Ito/Girsanov form of HMC (continuous time limit)                |
| midpt    | Finite time HMC using the mid point integrator                  |
| leapfrog | Finite time HMC using the LeapFrog (velocity Verlet) integrator |

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
