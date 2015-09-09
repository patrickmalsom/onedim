# HMC algorithm for a 1D potential

This code runs the HMC algorithm for a particle in a 1D potential.
The algorithm is implimented for the following methods:

|  method  | description |
| -------- | ----------- |
| ito      | Ito/Girsanov form of HMC (continuous time limit)                |
| midpt    | Finite time HMC using the mid point integrator                  |
| leapfrog | Finite time HMC using the LeapFrog (velocity Verlet) integrator |

## Running the simulation

`onedimHMC.py` is the main executable used to run the simulation.
Help using this program can be found by running `./onedimHMC.py -h`.

There are 2 required arguments:
1. **method**: integration method (ito, finite)
2. **HMC**: number of HMC loops to run

There are many optional arguments which may be changed as well, depending on the run to be performed, including

| argument | description |
| -------- | ----------- |
| --potential    | Potential def'n, found in `potential_defns` directory. See more **TODO** here. | 
| --inpath  | Input (starting) path. |
| --temp    | Configurational temperature. |
| --dt     | Time step along the path. |
| --dtau   | Time step between paths. |
| --pathlen    | Number of positions along the path. |
| --mdsteps     | Number of molecular dynamics steps to perform between each HMC step. |
| --writefiles | Number of intermediate paths to write to `output_paths` directory.|
| --seed | Seed (integer) for the random number generator. Default uses sysrandom. |
| --debug | Optional debugging file. defaults to off. See more **TODO** here. |

## Compiling the C library

This code uses python's ctypes to exploit the speed of a compiled C library.
Before running this code, the C library must be compiled using the `make` command.

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
