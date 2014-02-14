/*   GenPaths.c 
 *   Written Spring 2014 -- Patrick Malsom
 
program used to calculate the KL derivative stats.
takes m(t) and A(t) and generates expectation values.

This is just a shared library that is linked with guided.py
*/

// =========================================
// Definitions and Prototypes
// =========================================
//STD Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//OpenMP libraries
#include <omp.h>
//GNU Scientific Libraries
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

// Path Parameters
#define TEMP    0.15
#define NUMb    10001
#define DELTAt  0.001
#define XSTART  -0.97

//Stores m, A, and KL distance expectation values
typedef struct _averages
{
  double m;
  double A;
  double xbar;
  double xxbar;
  double expVal[5];
} averages;

// Function prototypes
void GenPaths(int loopIters, int RNGseed, averages* bead);
double genStep(double x0, double v0h, averages* bead, int n);
double energyChange(double x0, double x1, averages* bead, int n);
double pot(double x);
double DeltaU(double x, averages* bead, int n);


void GenPaths(int loopIters, int RNGseed, averages* bead)
{

  //arrays to store the random numbers in (generate a full path at once)
  double GaussRand[NUMb]; 
  double UniformRand[NUMb]; 

  int beadInc;
  double x0;
  double x1;

  int acc;
  int rej;
  int loop;

  double v0h;
  double v0hPref=sqrt(2.0*DELTAt*TEMP);
  double deltaE;

  double evXM;
  double evUU0;
  double evA;

  //===============================================================
  // GNU Scientific Library Random Number Setup
  //===============================================================
  const gsl_rng_type * RanNumType;
  gsl_rng *RanNumPointer;
  gsl_rng_env_setup();
  RanNumType = gsl_rng_default;
  RanNumPointer= gsl_rng_alloc (RanNumType);
  gsl_rng_set(RanNumPointer, RNGseed);
  int RNGcount;

  // =========================================
  // Main Loop
  // =========================================
  acc=0;
  rej=0;

  // Perform the sampling
  for(loop=0; loop<loopIters;loop++)
  {
    for(RNGcount=0;RNGcount<NUMb; RNGcount++){
      GaussRand[RNGcount]= gsl_ran_gaussian_ziggurat(RanNumPointer,1.0);
    }
    for(RNGcount=0;RNGcount<NUMb; RNGcount++){
      UniformRand[RNGcount]= gsl_rng_uniform(RanNumPointer);
    }

    x0=XSTART;
    x1=XSTART;
    
    // generating the single path
    for(beadInc=0; beadInc<NUMb; beadInc++)
    {
      x0=x1;
      v0h=v0hPref*GaussRand[beadInc];

      x1=genStep(x0, v0h, bead, beadInc);

      deltaE = energyChange(x0, x1, bead, beadInc);

      if(exp(-deltaE/TEMP) >= UniformRand[beadInc]){
        acc++;
      }
      else{
        rej++;
        x1=x0;
      }

      // accumulate expectation values and means
      evXM= x1-bead[beadInc+1].m; // x - m
      evUU0 = DeltaU(x1, bead, beadInc); // U - U0
      evA = bead[beadInc+1].A; // A

      bead[beadInc].xbar+=x1;
      bead[beadInc].xxbar+=x1*x1;
      bead[beadInc].expVal[0] += evUU0;
      bead[beadInc].expVal[1] += (-(evXM)*evA);
      bead[beadInc].expVal[2] += (-(evXM)*evA)*evUU0;
      bead[beadInc].expVal[3] += (0.5*evXM*evXM);
      bead[beadInc].expVal[4] += (0.5*evXM*evXM*evUU0);
    }
  }
}

// =========================================
// Functions
// =========================================

// =========================================
//genStep: 
double genStep(double x0, double v0h, averages* bead, int n)
{
  double ah0=bead[n].A*0.5;
  double ah1=bead[n+1].A*0.5;
  double m0=bead[n].m;
  double m1=bead[n+1].m;

  return m1 + (-m1+v0h+x0+ah0*(m0-x0)*DELTAt)/(1.0+ah1*DELTAt);
}

// =========================================
//energyChange: x0, x1, t0, t1, thalf
double energyChange(double x0, double x1, averages* bead, int n)
{
  double A0=bead[n].A;
  double A1=bead[n+1].A;
  double m0=bead[n].m;
  double m1=bead[n+1].m;

  return 0.5*( (-A0*(m0-x0)*(m0-x1)+A1*(m1-x0)*(m1-x1)) );

}

// =========================================
double pot(double x)
{
  // U = (x^2-1)^2
  //return 1.0 + x*x*(-2.0 + x*x);
  return (x*x - 1.0)*(x*x - 1.0);
}

// =========================================
//DeltaU: x m A
double DeltaU(double x, averages* bead, int n)
{
  double xminm= x - bead[n].m;
  return pot(x) - 0.5*bead[n].A*xminm*xminm;
}

// =========================================
