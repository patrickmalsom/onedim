#!/bin/bash
module load python/2.7 
module load libgfortran 
module load atlas 
module load lapack 
module load all-pkgs 
./forward_Bs_length.py "$@"
