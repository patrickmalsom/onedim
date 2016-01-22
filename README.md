# HMC algorithm for a 1D potential

This code runs the HMC algorithm for a particle in a 1D potential.

## Compiling the C libraries

#### Requirements
The code is compiled using `gcc` and openMP is a requirement to enable multi-threading.
The main python code is written in Python 2.7 (or 2.6) and requires numpy be installed for your version of python.

#### RHEL/CentOS

The yum package manager can handle all of the dependencies, simply install the following packages:

```
autotools
gcc
libgomp
python-numpy
```

#### Building

This code uses python's ctypes to exploit the speed of a compiled C library, which must be compiled before running.
Simply run make to compile the project:

```
make
```

## Running the simulation

`onedimHMC.py` is the main executable used to run the simulation.
Help using this program can be found by running `./onedimHMC.py -h`.

There are 2 required arguments:

1. **method**: integration method
    - **ito**: Ito/Girsanov form of HMC (continuous time limit)
    - **midpt**: Finite time HMC using the mid point integrator
    - **leapfrog**: Finite time HMC using the LeapFrog (velocity Verlet) integrator
    - **simpsons**: Finite time HMC using the LeapFrog (velocity Verlet) integrator
2. **HMC**: (integer) number of HMC loops to run

### Optional command line arguments

There are many optional arguments which may be changed depending on the run to be performed.
These optional arguments default to sane settings for the [Narrow-Broad potential](potential_defns/fatter_skinny.c).

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

### Example run

Running the simulation is performed using the python executable.
Lets try this out with the continuous-time integration method (*ito*) and a single HMC loop by running:

```
./onedimHMC.py ito 1
```

## Creating a new potential

Including another potential to the library is very simple.
The only requirement is that the potential be continuously differentiable (the force, F = -dV/dx, must be continuous).
The user must calculate and enter following functions into the new potential definition file:

|      function       | description | 
| ------------------- | ----------- |
| Pot(x)              | V(x): potential |
| Force(x)            | $F(x) = -\frac{d}{dx} V(x)$ |
| ForcePrime(x)       | $F'(x) = -\frac{d^2}{dx^2} V(x)$ |
| ForceDoublePrime(x) | $F''(x) = -\frac{d^3}{dx^3} V(x)$ |

These functions can be calculated using python and sympy (see the [broad_narrow potential example](potential_calculations/broad_narrow) ) with the following python code:

``` python
# sympy package is required
import sympy as sym

# define sympy symbols and functions
x=sym.symbols('x')
pot = sym.Function('pot')

# declare the potential function and print it in Horner Form
pot=((8 - 5*x)**8 * (2 + 5*x)**2)/2**26
sym.polys.polyfuncs.horner(sym.N(pot))

# calculate the force
F=sym.diff(-pot,x)
sym.polys.polyfuncs.horner(sym.N(F))

# calculate the force prime
Fp=sym.diff(-pot,x,2)
sym.polys.polyfuncs.horner(sym.N(Fp))

# calculate the force double prime
Fpp=sym.diff(-pot,x,3)
sym.polys.polyfuncs.horner(sym.N(Fpp))
```

Now, add these functions to a C file in the `potential_defns` directory.
This file may have any name, but must end with `.c` for the Makefile to be aware of it.

#### Required include files

The potential definition C file must include the following headers at the top of the file

``` c
#include "../onedimHMC_struct.h"
#include "../onedimHMC_forces.h"
#include "../onedimHMC.c"
```

Please refer to an existing potential definition file for an [example](potential_defns/fatter_skinny.c) of how to create a new potential.

