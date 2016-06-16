#!/usr/bin/env python

import numpy as np
import sys

data=np.loadtxt(sys.argv[1])
print "%s %i %f %f"%(sys.argv[1],len(data),np.average(data),np.std(data))

