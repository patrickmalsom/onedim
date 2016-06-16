#!/bin/bash

cat log/*.out > allout.txt

cat allout.txt | grep -v OVER | grep LeapFrog | sed s/LeapFrog//g| grep -v "0.0300" > alloutLeapFrog.txt
cat allout.txt | grep -v OVER | grep MidPt | sed s/MidPt//g| grep -v "0.0300" > alloutMidPt.txt
cat allout.txt | grep -v OVER | grep Simpson | sed s/Simpson//g| grep -v "0.0300" > alloutSimpson.txt
