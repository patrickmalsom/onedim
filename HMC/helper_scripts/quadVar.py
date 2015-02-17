#!/usr/bin/env python2.7

import numpy
import sys

# 1st arg: input path file
temp=numpy.loadtxt(sys.argv[1])
# 2nd arg: delta t step size
dt=float(sys.argv[2])

NumB=len(temp)

print "Conf Temp: %f" % (sum([ (temp[i] - temp[i+1])**2 for i in range(NumB-1)])/(2.*dt*NumB))
