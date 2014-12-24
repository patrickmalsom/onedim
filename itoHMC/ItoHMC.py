#!/usr/bin/env python

# Ito HMC

# ===============================================
# Import and setup libraries
# ===============================================
import numpy as np
import math
import sys
import random
import hashlib
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

# Argparse: command line options
import argparse
parser = argparse.ArgumentParser(description='Ito HMC algorithm (finite time step) for 1D external potential')
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
parser.add_argument('--deltatau', type=float, default=10**(-5), 
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
clib=ctypes.CDLL("ItoHMC.so")

# library functions
#c_Pot=clib.Pot
#clib.Pot.restype = DOUBLE
c_calcPotentials=clib.calcPotentials
c_calcLInverse=clib.LInverse

c_calcSPDEpos=clib.calcSPDEpos
c_calcMDpos=clib.calcMDpos
c_genBB=clib.generateBB

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
# IMPORTANT: These structs must EXACTLY MATCH the structs defined in the C code!

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
              ('bb', DOUBLE),
              ('F', DOUBLE),
              ('Fp', DOUBLE),
              ('Fpp', DOUBLE),
              ('G', DOUBLE),
              ('gradG', DOUBLE),
              ('LinvG', DOUBLE)
             ]
# define the array of averages that is to be passed to the C routine
pathType=Averages*NumB



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
    c_calcPotentials(pathOld, params);
    c_calcLInverse(pathOld, params);

def FillCurState():
    global pathCur, params
    c_calcPotentials(pathCur, params);
    c_calcLInverse(pathCur, params);

def FillNewState():
    global pathNew, params
    c_calcPotentials(pathNew, params);
    c_calcLInverse(pathNew, params);

def saveStartingState(pathSave):
    global pathCur
    for i in range(NumB):
        pathSave[i]=pathCur[i].pos


def printState(identifier):
    print "%s" % (identifier),
    print "\tqv: %1.6f" % c_quadVar(pathCur,params),
    print "\tDelta E: %1.15f" % (Echange*2.*eps) ,
    print "\tExp(-dE*dt/2/eps): %1.15f" % math.exp(-Echange)

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
    print '  input file name=%s' % args.infile
    print '  md5 hash: '+hashlib.md5(open(args.infile).read()).hexdigest()
    print '  eps  = %f' % eps
    print '  dt   = %e' % args.deltat
    print '  dtau = %e' % args.deltatau
    print '  Nb   = %d' % NumB
    print '  HMC  = %d' % args.HMC
    print '  MD   = %d' % args.MD
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
class TestClass:

    def test_setUp(self):
        print "starting the test suite"
        # dont run any HMC loops
        args.HMC=0
        assert 1

    def test_tearDown(self):
        assert 1
        pytest.exit("ending the tests")


# ===============================================
# MAIN loop
# ===============================================

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
    # generate the brownian bridge from the random numbers
    c_genBB(pathCur,params)

    ##====================================
    ##PinskiDebug
    #BBtemp=np.loadtxt("RNGBB.dat")
    #for i in range(10001):
    #    pathCur[i].bb=BBtemp[i]
    ##====================================

    # filled the entire path struct
    FillCurState()
    c_calcLInverse(pathCur,params)

    ##====================================
    ##PinskiDebug
    #np.savetxt("LinvG.dat",[pathCur[i].LinvG for i in range(10001)])
    ##====================================

    # calculate the new positions
    c_calcSPDEpos(pathCur, pathNew, params)

    # calculate all of the struct arrays
    FillNewState()
    c_calcLInverse(pathNew,params)

    # reset the energy change accumulator at the beginning of the SPDE step
    Echange=0.0
    Echange+=c_calcEnergyChange(pathCur,pathNew,params)

    rotatePaths()
    print "======== SPDE =========="
    printState("SPDE")
#    print "0: %f 1: %f mid: %f end: %f" % (pathCur[0].pos, pathCur[1].pos,pathCur[5000].pos, pathCur[NumB-1].pos)
#    print "BB0: %f 1: %f mid: %f end: %f" % (pathCur[0].bb, pathCur[1].bb,pathCur[5000].bb, pathCur[NumB-1].bb)
    print "========================"


    MDloops=max(1,int(args.MD*(0.5 + np.random.random())))

    ##====================================
    ##PinskiDebug
    #MDloops=int(args.MD)
    ##====================================

    for MDIter in range( MDloops ):
 
        # calculate the new positions
        c_calcMDpos(pathOld, pathCur, pathNew, params)


        # calculate all of the struct arrays
        FillNewState()
        c_calcLInverse(pathNew,params)

        Echange+=c_calcEnergyChange(pathCur,pathNew,params)

 
        # rotate the path structs
        rotatePaths()

        if MDIter % int(int(args.MD)/5.) == 0:
            printState("MDloop "+str(MDIter))

    printState("MDloop "+str(MDloops))

    # Metropolis Hasitings Step
    if math.exp(-Echange) > np.random.random():
    ##====================================
    ##PinskiDebug
    #if math.exp(-Echange) > 0.:
    ##====================================

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



    if HMCIter % max(1,int(int(args.HMC)/float(args.WriteFiles))) == 0:
        # write/plot the path and save to file
        writeCurPath("Path"+str(HMCIter)+".dat")
        #plt.plot([pathCur[i].pos for i in range(NumB)])
        #plt.savefig('testplot'+str(plotiter)+'.png')
        #plt.close()
        #makeHistogram([pathCur[i].pos for i in range(NumB)],str(plotiter))
        #plotiter+=1








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
