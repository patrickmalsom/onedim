#!/bin/bash

########################################################################
## Run on only on Open Science Grid
## Import gcc via cvmfs
#source /cvmfs/oasis.opensciencegrid.org/osg/modules/lmod/5.6.2/init/bash;
#module load gcc/4.6.2;
## Import python and numpy via cvmfs
#source /cvmfs/oasis.opensciencegrid.org/osg/modules/lmod/5.6.2/init/bash;
#module load python/2.7;
#module load atlas;
#module load lapack;
#module load all-pkgs;
########################################################################

# Export the current dir as a library path 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)

# run the python script and pass all of the options to argparse
./ItoHMC.py $@
