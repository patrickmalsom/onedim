#!/usr/bin/env python

# MALA routine

# ===============================================
# Import and setup libraries
# ===============================================
import numpy as np
import math
import sys
import random
import hashlib
import os
##plotting with text interface
#import matplotlib
#matplotlib.use('Agg')
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
    help='time step between paths;        default=10**(-5)')
parser.add_argument('--Num', type=int, default=10001, 
    help='path length (num beads);        default=10001')
parser.add_argument('--HMC', type=int, default=2, 
    help='HMC loops (SDE+MD+MC);          default=2')
parser.add_argument('--MD', type=int, default=150, 
    help='MD steps per full HMC;          default=150')
parser.add_argument('--WriteFiles', type=int, default=5, 
    help='Number of files to write;       default=5')
parser.add_argument('--RNGseed', type=int, 
    help='random number seed;             default=random')
args = parser.parse_args()

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
c_calcPosBar=clib.calcPosBar
c_calcForces=clib.calcForces
c_calcForcesBar=clib.calcForcesBar
c_calcForcesPrimeBar=clib.calcForcesPrimeBar
c_calcForcesDoublePrimeBar=clib.calcForcesDoublePrimeBar
c_calcDeltae=clib.calcDeltae
c_calcdg=clib.calcdg
c_calcSPDErhs=clib.calcSPDErhs
c_calcMDrhs=clib.calcMDrhs
c_calcPhi=clib.calcPhi

c_GaussElim=clib.GaussElim
# Gaussian elimination for computing L.x=b where L is the second deriv matrix

c_calcEnergyChange=clib.calcEnergyChange
clib.calcEnergyChange.restype = DOUBLE

c_quadVar=clib.quadVar
clib.quadVar.restype = DOUBLE


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
              ('posBar', DOUBLE),
              ('randlist', DOUBLE),
              ('F', DOUBLE),
              ('Fbar', DOUBLE),
              ('Fpbar', DOUBLE),
              ('Fppbar', DOUBLE),
              ('deltae', DOUBLE),
              ('dg', DOUBLE),
              ('Phi', DOUBLE),
              ('rhs', DOUBLE),
              ('bb', DOUBLE),
              ('Fp', DOUBLE),
              ('Fpp', DOUBLE),
              ('G', DOUBLE),
              ('gradG', DOUBLE),
              ('LinvG', DOUBLE)
             ]
# define the array of averages that is to be passed to the C routine
pathType=Averages*NumB

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
    c_calcPosBar(pathOld, params);
    c_calcForces(pathOld, params);
    c_calcForcesBar(pathOld, params);
    c_calcForcesPrimeBar(pathOld, params);
    c_calcForcesDoublePrimeBar(pathOld, params);
    c_calcDeltae(pathOld, params);
    c_calcdg(pathOld, params);
    c_calcPhi(pathOld, params);

def FillCurState():
    global pathCur, params
    c_calcPosBar(pathCur, params);
    c_calcForces(pathCur, params);
    c_calcForcesBar(pathCur, params);
    c_calcForcesPrimeBar(pathCur, params);
    c_calcForcesDoublePrimeBar(pathCur, params);
    c_calcDeltae(pathCur, params);
    c_calcdg(pathCur, params);
    c_calcPhi(pathCur, params);

def FillNewState():
    global pathNew, params
    c_calcPosBar(pathNew, params);
    c_calcForces(pathNew, params);
    c_calcForcesBar(pathNew, params);
    c_calcForcesPrimeBar(pathNew, params);
    c_calcForcesDoublePrimeBar(pathNew, params);
    c_calcDeltae(pathNew, params);
    c_calcdg(pathNew, params);
    c_calcPhi(pathNew, params);

def saveStartingState(pathSave):
    global pathCur
    for i in range(NumB):
        pathSave[i]=pathCur[i].pos


def printState(identifier):
    print "%s" % (identifier),
    print "qv: %1.6f" % c_quadVar(pathCur,params),
    print "\tDelta E: %1.6f" % Echange ,
    print "\tExp(-dE*dt/2/eps): %1.6f" % math.exp(-Echange*deltat/(2.0*eps))

def printTrans():
    global pathNew

    transCt=0
    xstart=pathNew[0].pos

    basinLeft=-2/3.
    basinRight=4/3.

    basin=basinLeft

    for i in xrange(NumB):
        if (pathNew[i].pos >4/3.) and (basin == basinLeft):
            basin=basinRight
        if (pathNew[i].pos < -2/3.) and (basin == basinRight):
            basin=basinLeft
            transCt+=1
    print "Transitions: %d" % (transCt)

def printParams():
    print '------------------------------------------------'
    print 'New HMC algorithm (finite time step) for 1D external potential'
    print '  input file : %s' % args.infile
    print '  md5 hash   : '+hashlib.md5(open(args.infile).read()).hexdigest()
    print '  quadvar eps: TODO'
    print '  eps  = %f' % eps
    print '  dt   = %e' % args.deltat
    print '  dtau = %e' % args.deltatau
    print '  Nb   = %d' % NumB
    print '  HMC  = %d' % args.HMC
    print '  MD   = %d' % args.MD
    print '  MD*h = %f' % (float(args.MD) * math.sqrt(2.0*args.deltatau))
    print "  seed = %d" % args.RNGseed
    print '------------------------------------------------'

def writeCurPath(fileName):
    global pathCur
    np.savetxt(fileName,np.array([pathCur[i].pos for i in range(NumB)]))

#def makeHistogram(path,pltname):

#    #calculate the partition function for use in the plots
#    expPot = lambda x: math.exp(-(((3.*x+2.)**2 * (3.*x-4.)**4)/1024.)/0.15)
#    Z0=1.1229622054
#    ### mathematica code to find partition function normalization
#    #NIntegrate[Exp[-((3 x + 2)^2 * (3 x - 4)^4/1024)/0.15], {x, -100, 100}]

#    HistData=[0 for i in range(400)]
#    for i in range(len(path)):
#        HistData[int((path[i]+1.5)*100)]+=1

#    invPathLen=1.0/(len(path)*0.005)
#    plt.plot([i for i in np.linspace(-1.5,2.5,400)],[expPot(i)/Z0 for i in np.linspace(-1.5,2.5,400)])
#    plt.plot([i for i in np.linspace(-1.5,2.5,400)],np.array(HistData)*invPathLen)
#    plt.savefig('Histplot'+pltname+'.png')
#    plt.close()

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

plotiter=100000

## start the profiler
#pr = cProfile.Profile()
#pr.enable()

printParams()

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
    printState("SPDE")
    print "========================"


    #MDloops=max(1,int(args.MD*(0.5 + np.random.random())))
    ##====================================
    ##PinskiDebug
    MDloops=int(args.MD)
    ##====================================

    #for MDIter in range( max(1,int(args.MD*(0.5 + np.random.random()))) ):
    for MDIter in range( MDloops ):
 
        # calculate the current state 
        c_calcMDrhs(pathOld, pathCur, params)

        # generate the pathNew positions
        c_GaussElim(pathCur,pathNew,params)

        # calculate all of the struct arrays
        FillNewState()

        Echange+=c_calcEnergyChange(pathCur,pathNew,params)

 
        # rotate the path structs
        rotatePaths()

        # print the MD state 
        # (less than 10 MD loops print all for debugging)
        if MDloops <= 10 and MDIter != MDloops-1:
            printState("MDloop "+str(MDIter))
        elif MDIter != MDloops-1 and MDIter % int(int(args.MD)/5.) == 0:
            printState("MDloop "+str(MDIter))

    printState("MDloop "+str(MDloops-1))


    # Metropolis Hasitings Step

    if math.exp(-Echange*deltat/(2.0*eps)) > np.random.random():
        # accept
        acc+=1
        posBasin=0
        for i in xrange(NumB):
            if pathCur[i].pos > 0:
                posBasin+=1
        print "posBasin: %i" % (posBasin)
        printTrans()
    else:
        rej+=1
        for i in range(NumB):
            pathCur[i].pos=savePath[i]

    print "acc: %d   rej: %d" % (acc,rej)



    if HMCIter % 10 == 0:
        # write/plot the path and save to file
        writeCurPath("testPath"+str(HMCIter)+".dat")
        #plt.plot([pathCur[i].pos for i in range(NumB)])
        #plt.savefig('testplot'+str(plotiter)+'.png')
        #plt.close()
        #makeHistogram([pathCur[i].pos for i in range(NumB)],str(plotiter))
        plotiter+=1








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
