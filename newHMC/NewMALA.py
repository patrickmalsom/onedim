#!/usr/bin/env python

# MALA routine

# ===============================================
# Import and setup libraries
# ===============================================
import numpy as np
import math
import sys
import random
#import matplotlib.pyplot as plt

# Argparse: command line options
import argparse
parser = argparse.ArgumentParser(description='New HMC algorithm (finite time step) for 1D external potential')
parser.add_argument('-t','--test', action="store_true",
    help='run py.test unit tests')
parser.add_argument('-i','--infile', type=str, default='inFile.dat', 
    help='input path positions;           default=inFile.dat')
parser.add_argument('-o','--outfile', type=str, default='outFile.dat', 
    help='output path positions;          default=outFile.dat')
parser.add_argument('-T','--temperature', type=float, default=0.15, 
    help='configurational temperature;    default=0.15')
parser.add_argument('--deltat', type=float, default=0.005, 
    help='time step along the path;       default=0.005')
parser.add_argument('--deltatau', type=float, default=10**(-7), 
    help='time step between paths;        default=10**(-7)')
parser.add_argument('--Num', type=int, default=10001, 
    help='path length (num beads);        default=10001')
parser.add_argument('--HMC', type=int, default=2, 
    help='HMC loops (SDE+MD+MC);          default=2')
parser.add_argument('--MD', type=int, default=150, 
    help='MD steps per full HMC;          default=100')
parser.add_argument('--RNGseed', type=int, 
    help='random number seed;             default=random')
args = parser.parse_args()

# numpy.random.seed(random.SystemRandom().randint(1,100000000000))


# python profiler
import cProfile, pstats, StringIO

# if testing mode flag is given
if args.test:
    # py.test unit testing suite
    import pytest
    pytest.main("NewMALA.py")

# ===============================================
# Ctypes
# ===============================================
import ctypes

# define some c variable types
DOUBLE=ctypes.c_double
PDOUBLE=ctypes.POINTER(DOUBLE)
INT=ctypes.c_int
PINT=ctypes.POINTER(INT)

# import the c code library
clib=ctypes.CDLL("NewMALA-func.so")

# library functions
#c_Pot=clib.Pot
#clib.Pot.restype = DOUBLE

c_calcForces=clib.calcForces
# fill out the Force variables for the struce

c_calcHessian=clib.calcHessian
# fill out the Hessian variables for the struce

c_calcdg=clib.calcdg
# calculate dG/dx for the structure

c_calcSPDErhs=clib.calcSPDErhs
c_calcMDrhs=clib.calcMDrhs
# calculate the right hand side vector for the SPDE


c_GaussElim=clib.GaussElim
# Gaussian elimination for computing L.x=b where L is the second deriv matrix

c_calcEnergyChange=clib.calcEnergyChange
clib.calcEnergyChange.restype = DOUBLE

c_quadVar=clib.quadVar
clib.quadVar.restype = DOUBLE

c_calcDeltae=clib.calcDeltae
c_calcPhi=clib.calcPhi


# ===============================================
# constants/parameters
# ===============================================
deltat=args.deltat # time step along the path
eps=args.temperature # configurational temperature
deltatau=args.deltatau # time step between paths
NumB=args.Num# number of 'beads' (total positions in path, path length)

invdt=1/deltat # reciprocal of time step
noisePref=math.sqrt(4.0*eps*deltatau*invdt) # scaling for the random noise term
r=deltatau*0.5*invdt*invdt # constant matrix element in L (I+-rL)
NumU=NumB-1 # number of time steps
NumL=NumB-2 # size of matrix (I+-rL)

# ===============================================
# C struct setup
# ===============================================

# save parameters to a struct to pass to the C routines
class Parameters(ctypes.Structure):
  _fields_ = [('deltat',DOUBLE),
              ('invdt', DOUBLE),
              ('eps', DOUBLE),
              ('deltatau', DOUBLE),
              ('noisePref', DOUBLE),
              ('r', DOUBLE),
              ('NumB', INT),
             ]
# define the array of averages that is to be passed to the C routine
paramType=Parameters

# make a class which is a clone of the C struct
# this structs pointer will be passed between C and Python code
class Averages(ctypes.Structure):
  _fields_ = [('pos', DOUBLE),
              ('randlist', DOUBLE),
              ('Force', DOUBLE),
              ('Hessian', DOUBLE),
              ('deltae', DOUBLE),
              ('dg', DOUBLE),
              ('Phi', DOUBLE),
              ('rhs', DOUBLE)
             ]
# define the array of averages that is to be passed to the C routine
pathType=Averages*NumB


# ===============================================
# Function declarations
# ===============================================
def Pot(x):
    return 1.+ x*x*(-3.375+x * (1.6875 +x * (2.84766 +(-2.84766+0.711914 * x) * x)))

def Force(x):
    return x * (6.75 + x * (-5.0625 + x * (-11.3906 + (14.2383 - 4.27148 * x) * x)))

def ForcePrime(x):
    return 6.75 + x * (-10.125 + x * (-34.1719 + (56.9531 - 21.3574 * x) * x))

def g(xm1,x0,x1):
    g1st=ForcePrime(x0)*Force((x0))
    g2nd=-0.5*invdt*(Force((x1))-2.0*Force((x0))+Force((xm1)))
    g3rd=-0.5*invdt*ForcePrime(x0)*(x1-2.0*x0+xm1)
    return g1st+g2nd+g3rd

def GaussElim(r,bVec):
    # Gaussian elimination for computing L.x=b where
    #   L is tridiagonal matrix
    #     mainDiag: 1+2r
    #     upper(lower)Diag: -r
    #   b is known vector: bVec is input
    #   x is unknown vector: returned at the end

    Num=len(bVec)
    
    al=[-r for i in range(Num-1)]
    am=[(1.0+2.0*r) for i in range(Num)]
    au=[-r for i in range(Num-1)]


    # Gaussian Elimination
    for i in range(Num-1):
        temp=-al[i]/am[i]
        #al[i]+=am[i]*temp
        am[i+1]+=au[i]*temp
        bVec[i+1]+=bVec[i]*temp

    # Back substitution
    for i in range(Num-1):
        temp=-au[Num-i-2]/am[Num-i-1]
        #au[Num-i-2]+=am[Num-i-1]*temp
        bVec[Num-i-2]+=bVec[Num-i-1]*temp
    

    # Divide by main diagonal
    for i in range(Num):
        bVec[i]=bVec[i]/am[i]
        #am[i]=am[i]/am[i]

    
    # Print the result
    return bVec


def rotatePaths():
    # rotate the structure pointers (old1->current0, current1->new0, new1->old0 )
    global pathOld,pathCur,pathNew
    temp=pathOld
    pathOld=pathCur
    pathCur=pathNew
    pathNew=temp
    del temp

def FillOldState():
    global pathOld, params
    c_calcForces(pathOld, params);
    c_calcHessian(pathOld, params);
    c_calcDeltae(pathOld, params);
    c_calcPhi(pathOld, params);
    c_calcdg(pathOld, params);

def FillCurState():
    global pathCur, params
    c_calcForces(pathCur, params);
    c_calcHessian(pathCur, params);
    c_calcDeltae(pathCur, params);
    c_calcPhi(pathCur, params);
    c_calcdg(pathCur, params);

def FillNewState():
    global pathNew, params
    c_calcForces(pathNew, params);
    c_calcHessian(pathNew, params);
    c_calcDeltae(pathNew, params);
    c_calcPhi(pathNew, params);
    c_calcdg(pathNew, params);

def saveStartingState(pathSave):
    global pathCur
    for i in range(NumB):
        pathSave[i]=pathCur[i].pos


def printState():
    print "qv: %1.6f" % c_quadVar(pathCur,params),
    print "\tDelta E: %1.6f" % Echange ,
    print "\tExp(-dE*dt/2/eps): %1.6f" % math.exp(-Echange*deltat/(2.0*eps))

# ===============================================
# Unit Tests
# ===============================================
def tFunc(x):
    return x+1

class TestClass:

    def test_setUp(self):
        print "starting the test suite"
        # dont run any HMC loops
        args.HMC=0
        assert 1

    def test_example_test1(self):
        assert tFunc(2)==3

    def test_example_test2(self):
        assert tFunc(2)==3 

    def test_tearDown(self):
        assert 1
        pytest.exit("ending the tests")

# Unit Tests
#def unitTest(path,Echange):
#    print "quadVar:" str(c_quadVar(path,params)) + "    Delta E:" + str(Echange)
    # print the quadratic variation
    #print str(c_quadVar(path,params))
    # verify that the boundary conditions have not changed
    #if inpath[0] != outpath[0] and inpath[-1] != outpath[-1]:
    #    sys.exit("Boundary Conditions Changed! Aborting...")
    
#def quadVar(poslist):
#    quadVarSum=sum([(poslist[i]-poslist[i+1])**2 for i in range(len(poslist)-1)])
#    quadVarSum/=(2.*eps*deltat*(len(poslist)-1))
#    print "quad var: " + str(quadVarSum)

## check for the 
#def checkBCs(path):
#    return [path[0],path[-1]]

#def printState(path0,path1,path2,beadNum):
#    print ""
#    print "================"
#    print "path0:"+str(path0[beadNum].pos)+"  "+str(path0[beadNum].Force)+"  "+str(path0[beadNum].dg)
#    print "path1:"+str(path1[beadNum].pos)+"  "+str(path1[beadNum].Force)+"  "+str(path1[beadNum].dg)
#    print "path2:"+str(path2[beadNum].pos)+"  "+str(path2[beadNum].Force)+"  "+str(path2[beadNum].dg)
#    print ""

#def plotPath(path):
#    tempPath=[path[i].pos for i in range(NumB)]
#    plt.plot(tempPath)
#    plt.show()

# ===============================================
# MAIN loop

if args.RNGseed is None:
  args.RNGseed = random.SystemRandom().randint(1,1000000)

np.random.seed(args.RNGseed)
print "RNG seed: %d" % args.RNGseed

inPath=np.loadtxt(args.infile)

# ===============================================
# Initializing the structs
#declare the params struct
params=paramType()
params.deltat=deltat
params.invdt=invdt
params.eps=eps
params.deltatau=deltatau
params.noisePref=noisePref
params.r=r
params.NumB=NumB

pathOld=pathType()
pathCur=pathType()
pathNew=pathType()


# output storage space
outPath=np.array([0.0 for i in np.arange(0,len(inPath),1)])
# save path (for rejection) storage space
savePath=np.array([0.0 for i in np.arange(0,len(inPath),1)])

# zero the acc/rej counters
acc=0
rej=0

## start the profiler
#pr = cProfile.Profile()
#pr.enable()

for i in np.arange(0,len(inPath),1):
    pathCur[i].pos=inPath[i]


for HMCIter in range(args.HMC):

    # save the current path to savePath in case of rejection
    for i in range(NumB):
        savePath[i]=pathCur[i].pos

    # noise vector
    for i in np.arange(1,len(inPath)-1,1):
        pathCur[i].randlist=np.random.normal(0,1)
    #printState(pathOld,pathCur,pathNew,1)
    randnumtemplist=[pathCur[i].randlist for i in range(NumB)]

    # calculate the current state 

    FillCurState()

    # filled the entire path struct
    c_calcSPDErhs(pathCur, params)

    # generates the pathCur positions
    c_GaussElim(pathCur,pathNew,params)

    # calculate all of the struct arrays
    FillNewState()

    # reset the energy change accumulator at the beginning of the SPDE step
    Echange=0.0
    Echange+=c_calcEnergyChange(pathCur,pathNew,params)


    #pTop=math.sqrt(params.deltatau/2.0) * ( (pathCur[i+1].pos-2.0*pathCur[i].pos+pathCur[i-1].pos) - pathCur[i].dg +(pathOld[i+1].pos-2.0*pathOld[i].pos+pathOld[i-1].pos) - pathOld[i].dg ) + params.noisePref*pathOld[i].randlist

             




    rotatePaths()
    print "======== SPDE =========="
    printState()
    print "========================"


    for MDIter in range(args.MD):
 
        # calculate the current state 
        c_calcMDrhs(pathOld, pathCur, params)

        # generate the pathNew positions
        c_GaussElim(pathCur,pathNew,params)

        # calculate all of the struct arrays
        FillNewState()

        Echange=c_calcEnergyChange(pathCur,pathNew,params)

 
        # rotate the path structs
        rotatePaths()

        if MDIter % 10 == 0:
            printState()


    # Metropolis Hasitings Step

    if math.exp(-Echange*deltat/(2.0*eps)) > np.random.random():
        # accept
        acc+=1
    else:
        rej+=1
        for i in range(NumB):
            pathCur[i].pos=savePath[i]

    print "acc: %d   rej: %d" % (acc,rej)








#dgTemp=[pathOld[i].dg for i in range(NumB)]
#plt.plot(dgTemp)
#plt.show()

#hessTemp=[pathOld[i].Hessian for i in range(NumB)]
#plt.plot(hessTemp)
#plt.show()

## profiler stop
#pr.disable()
#s = StringIO.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
#print s.getvalue()

# test to see if the path is still correct
for i in np.arange(0,len(inPath),1):
    outPath[i]=pathCur[i].pos

#if (np.loadtxt("outPathfive.dat") == outPath).all():
#    print "SUCCESS!"
#else:
#    print "FAIL"

np.savetxt(args.outfile,outPath)
