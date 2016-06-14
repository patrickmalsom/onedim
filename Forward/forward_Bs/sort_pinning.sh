#!/bin/bash

#cat log/*.out | grep startpin | sed 's/startpin\ Nb//g' > startpin_Bs.out
#cat log/*.out | grep startrevpin | sed 's/startrevpin\ Nb//g' > startrevpin_Bs.out
#cat log/*.out | grep bothpin | sed 's/bothpin\ Nb//g' > bothpin_Bs.out
#cat log/*.out | grep bothrevpin | sed 's/bothrevpin\ Nb//g' > bothrevpin_Bs.out

function printstats {
  cat startpin_Bs.out | grep "$1 " | cut -f2 -d " " > startpin_Bs_$1.out
  cat startrevpin_Bs.out | grep "$1 " | cut -f2 -d " " > startrevpin_Bs_$1.out
  cat bothpin_Bs.out  | grep "$1 " | cut -f2 -d " " > bothpin_Bs_$1.out
  cat bothrevpin_Bs.out  | grep "$1 " | cut -f2 -d " " > bothrevpin_Bs_$1.out
  ./calc_moments.py  startpin_Bs_$1.out
  ./calc_moments.py  startrevpin_Bs_$1.out
  ./calc_moments.py  bothpin_Bs_$1.out
  ./calc_moments.py  bothrevpin_Bs_$1.out
}

printstats 1000
printstats 2000
printstats 4000
printstats 6000
printstats 8000
printstats 10000
printstats 15000
printstats 20000
printstats 30000
printstats 40000
printstats 60000
printstats 80000
printstats 100000
