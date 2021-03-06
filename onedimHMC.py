#!/usr/bin/env python

# HMC 1D routine

##################### OVERVIEW #######################
# Ito                common             Finite
# ---                ------             ------
#                    Import modules
#                    argparse
#                    ctypes
#                    constant defns
#                    struct defns
#                    functions
#                      - rorate paths
#                      - fill states (x3)
#                      - saveStartState
#                      - print funcs
#                      - writeCurPath
#       ************* MAIN LOOP **************
#                    setRNGseed
#                    readInPath
#                    declare consts
#                    declare structs
#                    print params
#                    pathCur.pos=inPath
#       ************* HMC LOOP **************
#                    save path
#                <-- gen noise -->
# gen BB   
# fillCurIto                            fillCurFinite
# calcSPDEItoPos                        calcSPDEFinitePos
# ---                                       - calcSPDEFiniteRHS
# ---                                       - Gaussian Elim
# FillNewIto                            fillNewFinite
# EchangeIto -->                    <-- EchangeFinite
#                    rotateStruct
#                    printState
#       ************* MD LOOP ***************
#                <-- random MD loops -->
# calcMDItoPos                          calcMDFinitePos
# ---                                       - calcMDFiniteRHS
# ---                                       - Gaussian Elim
# fillNewIto                            FillNewFinite
# EchangeIto -->                    <-- EchangeFinite
#                    rotateStruct
#                    printState
#       ************* MHMC  ***************
#                    if(-dE>rand):
#                      acc
#                    else:
#                      rej
#                    writePaths
#       ************* End MAIN LOOP ***************
#                    write final path

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

parser = argparse.ArgumentParser(
        # argparse options
        description='Path-Space HMC algorithm for 1D potential',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

# Required arguments
parser.add_argument('method', type=str,
    help='HMC method to use;     options:{ito,midpt,leapfrog,simpsons}')
parser.add_argument('HMCsteps', type=int, default=2, 
    help='HMC loops (SDE+MD+MC)')

# Optional arguments
parser.add_argument('--potential', type=str, default='fatter_skinny',
    help='Potential to use;   Found in potential_defns dir')
parser.add_argument('--inpath', type=str, default='input_paths/fatterSkinny-T0p25-dt0p005-Nb30k-healed.dat', 
    help='input path positions')
parser.add_argument('--temp', type=float, default=0.25, 
    help='configurational temperature')
parser.add_argument('--dt', type=float, default=0.005, 
    help='time step along the path')
parser.add_argument('--dtau', type=float, default=10**(-7), 
    help='time step between paths')
parser.add_argument('--pathlen', type=int, default=30001, 
    help='path length (num beads)')
parser.add_argument('--mdsteps', type=int, default=500, 
    help='MD steps per full HMC')
parser.add_argument('--writefiles', type=int, default=0, 
    help='Number of files to write')
parser.add_argument('--seed', type=int, 
    help='random number seed; default uses SysRandom')
parser.add_argument('--debug', type=str, default='0',
    help='debug: save struct to file;     enter name of debug file to turn on')

#sys.argv includes a list of elements starting with the program
if len(sys.argv) < 2:
    parser.print_usage()
    sys.exit(0)

args = parser.parse_args()


# python profiler
import cProfile, pstats, StringIO

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
if os.path.isfile("potential_defns/"+args.potential+".so") is False:
    print "ERROR: Potential file does not exist! Exiting..."
    sys.exit(1)

clib=ctypes.CDLL("potential_defns/"+args.potential+".so")


# fill average position 
c_calcPosBar=clib.calcPosBar

# fill Forces
c_calcPot=clib.calcPot
c_calcForces=clib.calcForces
c_calcForcesBar=clib.calcForcesBar
c_calcForcesPrime=clib.calcForcesPrime
c_calcForcesPrimeBar=clib.calcForcesPrimeBar
c_calcForcesDoublePrime=clib.calcForcesDoublePrime
c_calcForcesDoublePrimeBar=clib.calcForcesDoublePrimeBar

# fill the Ito path potential
c_calcG=clib.calcG
c_calcgradG=clib.calcgradG

# Finite library functions
c_calcDeltae=clib.calcDeltae
c_calcdg=clib.calcdg
c_calcSPDErhs=clib.calcSPDEFiniterhs
c_calcMDrhs=clib.calcMDFiniterhs
c_calcPhi=clib.calcPhi
# Gaussian elimination for computing L.x=b where L is the second deriv matrix
c_GaussElim=clib.GaussElim

# Ito library functions
c_calcLInverse=clib.LInverse
c_calcSPDEItopos=clib.calcSPDEItopos
c_calcMDItopos=clib.calcMDItopos
c_genBB=clib.generateBB

# Finite energy chage calculation (returns double)
c_EChangeFinite=clib.calcEnergyChangeFinite
clib.calcEnergyChangeFinite.restype = DOUBLE
# Ito energy chage calculation (returns double)
c_calcEChangeIto=clib.calcEChangeIto
clib.calcEChangeIto.restype = DOUBLE

# quadratic variation of a path
c_quadVar=clib.quadVar
clib.quadVar.restype = DOUBLE


# ===============================================
# constants/parameters
# ===============================================
deltat=args.dt # time step along the path
eps=args.temp # configurational temperature
deltatau=args.dtau # time step between paths
NumB=args.pathlen# number of 'beads' (total positions in path, path length)

# finite method initialization
#   this will store the integer associated with the method to the 
#   parameters struct for use in the C code (deltae, Phi, dg)
availableMethods={'ito':0,'midpt':1,'leapfrog':2,'simpsons':3}
if args.method in availableMethods:
    methodInt=availableMethods[args.method]
else:
    print "Method does not exist! Choices are {ito,midpt,leapfrog,simpsons}. Exiting..."
    sys.exit(1)

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
              ('method', INT)
             ]
# define the array of averages that is to be passed to the C routine
paramType=Parameters

# make a class which is a clone of the C struct
# this structs pointer will be passed between C and Python code
class Averages(ctypes.Structure):
  _fields_ = [('pos', DOUBLE),
              ('posBar', DOUBLE),
              # random gaussian numbers
              ('randlist', DOUBLE),
              # forces
              ('U', DOUBLE),
              ('F', DOUBLE),
              ('Fbar', DOUBLE),
              ('Fp', DOUBLE),
              ('Fpbar', DOUBLE),
              ('Fpp', DOUBLE),
              ('Fppbar', DOUBLE),
              # Finite
              ('deltae', DOUBLE),
              ('dg', DOUBLE),
              ('Phi', DOUBLE),
              ('rhs', DOUBLE),
              # Ito
              ('bb', DOUBLE),
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

def FillCurIto():
    global pathCur, params
    c_calcPot(pathCur, params);
    c_calcForces(pathCur, params);
    c_calcForcesPrime(pathCur, params);
    c_calcForcesDoublePrime(pathCur, params);
    c_calcG(pathCur,params);
    c_calcgradG(pathCur,params);
    c_calcLInverse(pathCur, params);

def FillNewIto():
    global pathNew, params
    c_calcPot(pathNew, params);
    c_calcForces(pathNew, params);
    c_calcForcesPrime(pathNew, params);
    c_calcForcesDoublePrime(pathNew, params);
    c_calcG(pathNew,params);
    c_calcgradG(pathNew,params);
    c_calcLInverse(pathNew, params);

def FillCurFinite():
    global pathCur, params
    c_calcPosBar(pathCur, params);
    c_calcPot(pathCur, params);
    c_calcForces(pathCur, params);
    c_calcForcesBar(pathCur, params);
    c_calcForcesPrime(pathCur, params);
    c_calcForcesPrimeBar(pathCur, params);
    c_calcForcesDoublePrime(pathCur, params);
    c_calcForcesDoublePrimeBar(pathCur, params);
    c_calcDeltae(pathCur, params);
    c_calcdg(pathCur, params);
    c_calcPhi(pathCur, params);

def FillNewFinite():
    global pathNew, params
    c_calcPosBar(pathNew, params);
    c_calcPot(pathNew, params);
    c_calcForces(pathNew, params);
    c_calcForcesBar(pathNew, params);
    c_calcForcesPrime(pathNew, params);
    c_calcForcesPrimeBar(pathNew, params);
    c_calcForcesDoublePrime(pathNew, params);
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
    print "\tqv: %1.6f" % c_quadVar(pathCur,params),
    print "\tDelta E: %1.15f" % (Echange*2.*eps/deltat) ,
    print "\tExp(-dE*dt/2/eps): %1.15f" % math.exp(-Echange)

#def printTrans():
#    global pathNew

#    transCt=0
#    xstart=pathNew[0].pos

#    basinLeft=-2/3.
#    basinRight=4/3.

#    basin=basinLeft

#    for i in xrange(NumB):
#        if (pathNew[i].pos >4/3.) and (basin == basinLeft):
#            basin=basinRight
#        if (pathNew[i].pos < -2/3.) and (basin == basinRight):
#            basin=basinLeft
#            transCt+=1
#    print "Transitions: %d" % (transCt)

def printParams(path,params):
    print '--------------------------------------------------'
    print 'Path Space HMC algorithm for 1D external potential'
    print '--------------------------------------------------'
    print '  method : %s' % args.method
    print '  input file : %s' % args.inpath
    print '  md5 hash   : '+hashlib.md5(open(args.inpath).read()).hexdigest()
    print '  quadvar eps: %f' % ( quadraticVariation(path,params) )
    print '  eps  = %f' % params.eps
    print '  dt   = %e' % params.deltat
    print '  dtau = %e' % params.deltatau
    print '  Nb   = %d' % params.NumB
    print '  HMC  = %d' % args.HMCsteps
    print '  MD   = %d' % args.mdsteps
    print '  MD*h = %f' % (float(args.mdsteps) * math.sqrt(2.0*params.deltatau))
    print "  seed = %d" % args.seed
    print '------------------------------------------------'

def quadraticVariation(path,params):
    sumxx = sum([ (path[i+1].pos-path[i].pos)**2 for i in range(params.NumB-1) ])
    return sumxx/2.0/params.deltat/params.NumB

def writeCurPath(fileName):
    global pathCur
    np.savetxt("output_paths/"+fileName,np.array([pathCur[i].pos for i in range(NumB)]))

def initializeParams(params):
    params.deltat=deltat
    params.invdt=invdt
    params.eps=eps
    params.deltatau=deltatau
    params.noisePref=noisePref
    params.r=r
    params.NumB=NumB
    params.method=methodInt

def setRNGseed():
  # if there is no rng set on cmd line, generate a random one
  if args.seed is None:
    args.seed = random.SystemRandom().randint(1,1000000)
  # set the random seed of numpy
  np.random.seed(args.seed)

def calcSPDEFinitePos(path0,path1,params):
    # calculate the rhs vecor in the notes
    # the struct MUST be filled before this step
    c_calcSPDErhs(path0, params)

    # generates the pathCur positions by doing 
    # gaussian elim on the rhs vector
    c_GaussElim(path0,path1,params)

def calcMDFinitePos(path0,path1,path2,params):
    # calculate the current state 
    c_calcMDrhs(path0, path1, params)

    # generate the pathNew positions
    c_GaussElim(path1,path2,params)

def printStateMD(MDloops):
    # (less than 10 MD loops print all for debugging)
    if MDloops <= 10 and MDIter != MDloops-1:
        printState("MDloop "+str(MDIter))
    elif MDIter != MDloops-1 and MDIter % int(int(args.mdsteps)/5.) == 0:
        printState("MDloop "+str(MDIter))

def printDebugFinite(path0,path1,path2,params):
    # debugging function. prints many random values
        print "path0[20000]: de:%f Phi:%f dg:%f" % (path0[20000].deltae, path0[20000].Phi, path0[20000].dg)
        print "path1[20000]: de:%f Phi:%f dg:%f" % (path1[20000].deltae, path1[20000].Phi, path1[20000].dg)
        print "path2[20000]: de:%f Phi:%f dg:%f" % (path2[20000].deltae, path2[20000].Phi, path2[20000].dg)
        print "total deltae: %f" % (float(sum([path1[i].deltae for i in range(NumB)])))
        print "total abs(deltae): %f" % (float(sum([abs(path1[i].deltae) for i in range(NumB)])))
        print "total Phi: %f" % (float(sum([path1[i].Phi for i in range(NumB)])))
        print "total abs(Phi): %f" % (float(sum([abs(path1[i].Phi) for i in range(NumB)])))
        print "F1:%f F0:%f x1:%f x0:%f" % (path1[10001].F,path1[20000].F,path1[10001].pos,path1[20000].pos)

def saveDebugFinite(fileName,path0,path1,path2,params):
    np.savetxt(fileName+'.out', 
        [[
          pathCur[i].pos,
          pathCur[i].posBar,
          pathCur[i].randlist,
          pathCur[i].U,
          pathCur[i].F,
          pathCur[i].Fbar,
          pathCur[i].Fp,
          pathCur[i].Fpbar,
          pathCur[i].Fpp,
          pathCur[i].Fppbar,
          pathCur[i].deltae,
          pathCur[i].dg,
          pathCur[i].Phi,
          pathCur[i].rhs,
          pathCur[i].bb,
          pathCur[i].G,
          pathCur[i].gradG,
          pathCur[i].LinvG
         ] for i in range(NumB)]
        ,fmt='%8.8f')

def printDebugIto(path0,path1,path2,params):
    print "path0[20000]: de:%f G:%f gradG:%f" % (path0[20000].deltae, path0[20000].G, path0[20000].gradG)
    print "path1[20000]: de:%f G:%f gradG:%f" % (path1[20000].deltae, path1[20000].G, path1[20000].gradG)
    print "path2[20000]: de:%f G:%f gradG:%f" % (path2[20000].deltae, path2[20000].G, path2[20000].gradG)
    print "total deltae: %f" % (float(sum([path1[i].deltae for i in range(NumB)])))
    print "total abs(deltae): %f" % (float(sum([abs(path1[i].deltae) for i in range(NumB)])))
    print "total G: %f" % (float(sum([path1[i].G for i in range(NumB)])))
    print "total abs(G): %f" % (float(sum([abs(path1[i].G) for i in range(NumB)])))

def MHMC_test(Echange):
    global pathCur
    global pathNew
    global acc
    global rej
    if math.exp(-Echange) > np.random.random():
        # accept
        acc+=1
    else:
        # reject
        rej+=1
        # set current path postions to the saved path
        for i in range(NumB):
            pathCur[i].pos=savePath[i]
    print "acc: %d   rej: %d" % (acc,rej)
    
def printPosBasin():
    posBasin=0
    for i in xrange(NumB):
        if pathCur[i].pos > 0:
            posBasin+=1
    print "posBasin: %i" % (posBasin)

# ===============================================
# MAIN loop
# ===============================================
# set the seed for the random number generator
setRNGseed()

# read the input path specified on the cmd line
inPath=np.loadtxt(args.inpath)

# Initialize the structs for params and paths
#   need one parameter struct (constants)
#   need 3 path structs (old cur new)
params=paramType()
pathOld=pathType()
pathCur=pathType()
pathNew=pathType()


# fill the params struct with the run parameters
initializeParams(params)

# set the positions in pathCur to be inPath positions
for i in np.arange(0,len(inPath),1):
    pathCur[i].pos=inPath[i]

# print the run parameters
printParams(pathCur,params)

# output storage space
outPath=np.array([0.0 for i in np.arange(0,len(inPath),1)])
# save path (for rejection) storage space
savePath=np.array([0.0 for i in np.arange(0,len(inPath),1)])

# zero the acc/rej counters
acc=0
rej=0

# ============ SPDE/HMC LOOP =================
for HMCIter in range(args.HMCsteps):


    # save the current path to savePath in case of rejection
    for i in range(NumB):
        savePath[i]=pathCur[i].pos

    # generate the noise vector (random gaussian distributed list)
    for i in np.arange(1,len(inPath)-1,1):
        pathCur[i].randlist=np.random.normal(0,1)

    # ITO SPDE LOOP
    if (args.method == "ito"):
        # generate the brownian bridge from the random numbers
        c_genBB(pathCur,params)

        # filled the entire path struct
        FillCurIto()

        # calculate the new positions
        c_calcSPDEItopos(pathCur, pathNew, params)

        # calculate all of the struct arrays
        FillNewIto()

        # reset the energy change accumulator at the beginning of the SPDE step
        Echange=0.0
        Echange+=c_calcEChangeIto(pathCur,pathNew,params)

    # Finite SPDE LOOP
    else:
        # calculate the current state 
        FillCurFinite()

        # calculate the new position vector using the current state
        calcSPDEFinitePos(pathCur,pathNew,params)

        # calculate all of the struct arrays
        FillNewFinite()

        # reset the energy change accumulator at the beginning of the SPDE step
        Echange=0.0
        Echange+=c_EChangeFinite(pathCur,pathNew,params)

    rotatePaths()
    printState("SPDE")

    # ============ MD LOOP =================
    MDloops=max(1,int(args.mdsteps*(0.5 + np.random.random())))
    #MDloops=int(args.mdsteps)

    # ITO MD LOOP
    if (args.method == "ito"):

        for MDIter in range( MDloops ):

            # calculate the new positions
            c_calcMDItopos(pathOld, pathCur, pathNew, params)

            # calculate all of the struct arrays
            FillNewIto()

            Echange+=c_calcEChangeIto(pathCur,pathNew,params)

            # rotate the path structs
            rotatePaths()

            # print some of the MD states (~10 total)
            printStateMD(MDloops)
            if args.debug != '0':
                printDebugIto(pathOld,pathCur,pathNew,params)

    # Finite MD LOOP
    else:
        for MDIter in range( MDloops ):
            # calculate the new postions in pathNew
            calcMDFinitePos(pathOld,pathCur,pathNew,params)

            # calculate all of the struct arrays in pathNew
            FillNewFinite()

            Echange+=c_EChangeFinite(pathCur,pathNew,params)

            # rotate the path structs
            rotatePaths()

            # print some of the MD states (~10 total)
            printStateMD(MDloops)
            if args.debug != '0':
                printDebugFinite(pathOld,pathCur,pathNew,params)

    # print the run state for the final MD step
    printState("MDloop "+str(MDloops-1))

    # Metropolis Hasitings Monte-Carlo
    MHMC_test(Echange)
    printPosBasin()


    if args.writefiles is not 0:
        if (HMCIter % max(1,int(int(args.HMCsteps)/float(args.writefiles))) == 0) or (HMCIter == args.HMCsteps):
            writeCurPath("outPath"+str(HMCIter)+".dat")
    elif (HMCIter == args.HMCsteps):
        writeCurPath("outPath"+str(HMCIter)+".dat")

    sys.stdout.flush()

# struct debugging
if args.debug != '0':
    saveDebugFinite(args.debug,pathOld,pathCur,pathNew,params)
