#!/bin/bash

cat log/*.out | sed 's/LeapFrog/1/g' | sed 's/MidPt/2/g' | sed 's/Simpson/3/g' | grep -v -i overflow > stats.out
