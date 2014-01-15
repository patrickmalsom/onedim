#!/bin/bash

echo "Setting up Anaconda from OASIS..."
source /cvmfs/oasis.opensciencegrid.org/osg/palms/setup.sh
palmsdosetup anaconda

./SDE-statistics.py "$@"
