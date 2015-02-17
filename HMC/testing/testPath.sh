#!/usr/bin/env bash

make clean 
make

shopt -s expand_aliases
source ~/.bashrc

./onedimHMC.py --method ito --potential fatter_skinny -i input_paths/fatterSkinny-T0p25-dt0p005-Nb30k-end-healed.dat -T 0.25 --deltat 0.005 --deltatau 0.00000001 --Num 30001 --HMC 2 --MD 100 --RNGseed 10
awk '{getline t<"testing/outFileTestingIto.dat"; print $0-t}' output_paths/outPathFinal.dat > testing/outFileDiff.dat
cat testing/outFileDiff.dat | termplot

./onedimHMC.py --method finite --potential fatter_skinny -i input_paths/fatterSkinny-T0p25-dt0p005-Nb30k-end-healed.dat -T 0.25 --deltat 0.005 --deltatau 0.000001 --Num 30001 --HMC 2 --MD 100 --RNGseed 10
awk '{getline t<"testing/outFileTestingFinite.dat"; print $0-t}' output_paths/outPathFinal.dat > testing/outFileDiff.dat
cat testing/outFileDiff.dat | termplot
