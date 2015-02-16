#!/usr/bin/env bash

make clean 
make

shopt -s expand_aliases
source ~/.bashrc

./ItoHMC.py -i fatterSkinny-T0p25-dt0p005-Nb30k-end-healed.dat -T 0.25 --deltat 0.005 --deltatau 0.00000001 --Num 30001 --HMC 2 --MD 100 --RNGseed 10

awk '{getline t<"outFileTesting.dat"; print $0-t}' outPathFinal.dat > outFileDiff.dat
cat outFileDiff.dat | termplot
