#!/usr/bin/env bash

DT=0.005
EPS=0.25
XSTART=-0.4
NUM=100
METHOD=1 # {0:leapfrog, 1:midpt, 2:simpson}
MHMC=0
LOOPS=1
PRINT=1 # print run parameters on command line

usage() { echo "$0 usage:     [-t dt] [-e eps] [-x xstart] [-n num] [-m method] [-c MHMC] [-l loops]" && echo "" && echo "Detailed help:" && grep " .)\ #" $0 | sed s/\ \#//g ; exit 1; }

#usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }

while getopts "t:e:x:n:m:c:l:p" arg; do
    case "${arg}" in
        t) # time step - dt : default 0.005
            DT=${OPTARG} ;;
        e) # Temperature - eps : default 0.25
            EPS=${OPTARG} ;;
        x) # Initial position - x(0) : default -0.4
            XSTART=${OPTARG} ;;
        n) # Number of steps - num : default 100
            NUM=${OPTARG} ;;
        m) # quadrature {0:leapfrog,1:midpt,2:simpson} : default 1
            METHOD=${OPTARG} ;;
        c) # MHMC boolean {0:off,1:on} : default 0
            MHMC=${OPTARG} ;;
        l) # MHMC boolean {0:off,1:on} : default 0
            LOOPS=${OPTARG} ;;
        p) # suppress run information to stdout
            PRINT=0 ;;
        *)
            usage ;;
    esac
done
shift $((OPTIND-1))

GSL_RNG_SEED=$(shuf -i 2-999999999 -n 1) ./calc_move.out $DT $EPS $XSTART $NUM $METHOD $MHMC $LOOPS $PRINT 2> /dev/null
