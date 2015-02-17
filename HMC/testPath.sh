#!/usr/bin/env bash

make clean 
make

shopt -s expand_aliases
source ~/.bashrc

./onedimHMC.py --method ito -i fatterSkinny-T0p25-dt0p005-Nb30k-end-healed.dat -T 0.25 --deltat 0.005 --deltatau 0.00000001 --Num 30001 --HMC 2 --MD 100 --RNGseed 10
awk '{getline t<"outFileTestingIto.dat"; print $0-t}' outPathFinal.dat > outFileDiff.dat
cat outFileDiff.dat | termplot

./onedimHMC.py --method finite -i fatterSkinny-T0p25-dt0p005-Nb30k-end-healed.dat -T 0.25 --deltat 0.005 --deltatau 0.000001 --Num 30001 --HMC 2 --MD 100 --RNGseed 10
awk '{getline t<"outFileTestingFinite.dat"; print $0-t}' outPathFinal.dat > outFileDiff.dat
cat outFileDiff.dat | termplot
