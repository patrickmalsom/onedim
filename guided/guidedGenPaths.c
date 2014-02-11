/*   guidedGenPaths.c 
 *   Written Spring 2014 -- Patrick Malsom
 
 *  Outline of entire C routine: should be as general as possible
cmd line: ./a.out numLoops randSeed outfileName mean.dat A.dat
returns: file named "outfileName-numLoops.dat" that contains all stats

1. read in m and A
2. set up random number generator (read about how to do the seeds)
3. perform the sampling while keeping track of statistics
4. write out a file with numLoops in file name that have the statistics
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
//#define XSTART  -1.0


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
double genStep(double x0, double v0h, averages* bead, int n);
double energyChange(double x0, double x1, averages* bead, int n);
double pot(double x);
double DeltaU(double x, averages* bead, int n);
void writeConfig(averages* bead, int argc, char* argv[]);


// BEGIN MAIN FUNCTION
int main(int argc, char *argv[])
{

  averages *bead = calloc(NUMb,sizeof(averages));
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
  //printf("v0hpref:%+0.15e\n",v0hPref);
  double deltaE;

  double evXM;
  double evUU0;
  double evA;

  // =========================================
  // Reading in the input control paths (m and A)
  //==========================================
  int lineNum;

  lineNum = 0;
  FILE *fptrm = fopen(argv[4],"r");
  while( EOF != fscanf(fptrm,"%lf", &(bead[lineNum].m)) ) 
  {
    lineNum++;
  }

  lineNum = 0;
  FILE *fptrA = fopen(argv[5],"r");
  while( EOF != fscanf(fptrA,"%lf", &(bead[lineNum].A)) ) 
  {
    lineNum++;
  }

  //===============================================================
  // GNU Scientific Library Random Number Setup
  //===============================================================
  // Example shell command$ GSL_RNG_SEED=123 ./a.out
  const gsl_rng_type * RanNumType;
  gsl_rng *RanNumPointer;
  gsl_rng_env_setup();
  RanNumType = gsl_rng_default;
  RanNumPointer= gsl_rng_alloc (RanNumType);
  int RNGcount;
  //printf("Random Number Generator Type: %s \n", gsl_rng_name(RanNumPointer));
  //printf("RNG Seed: %li \n", gsl_rng_default_seed);

  // =========================================
  // Main Loop
  // =========================================


  acc=0;
  rej=0;

  // Perform the sampling
  for(loop=0; loop<atoi(argv[1]);loop++)
  {
    for(RNGcount=0;RNGcount<NUMb; RNGcount++){
      GaussRand[RNGcount]= gsl_ran_gaussian_ziggurat(RanNumPointer,1.0);
    }
    for(RNGcount=0;RNGcount<NUMb; RNGcount++){
      UniformRand[RNGcount]= gsl_rng_uniform(RanNumPointer);
    }

    //x0=XSTART;
    //x1=XSTART;
    x0=-0.9;
    x1=-0.9;
    
    // generating the single path
    for(beadInc=0; beadInc<NUMb; beadInc++)
    {
      x0=x1;
      v0h=v0hPref*GaussRand[beadInc];
      //printf("x1 : %+0.10e    v0h: %0.10e \n",x1,v0h);
      //printf("x1 : %0.10e     m: %0.10e \n",x1,bead[beadInc].m);

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
 
  //printf("acc:%d   rej:%d \n",acc,rej);
  // write out the file
  writeConfig(bead, argc, argv);

  return(0);
}

// =========================================
// Functions
// =========================================

// =========================================
//genStep: 
double genStep(double x0, double v0h, averages* bead, int n)
{
//  return x0 + v0h + DELTAt * ( bead[n+1].A*bead[n+1].m +bead[n].A*(bead[n].m - x0)) / (1.0 + DELTAt * bead[n+1].A);
  double ah0=bead[n].A*0.5;
  double ah1=bead[n+1].A*0.5;
  double m0=bead[n].m;
  double m1=bead[n+1].m;

  return m1 + (-m1+v0h+x0+ah0*(m0-x0)*DELTAt)/(1.0+ah1*DELTAt);
  //return ( x0 + v0h + 0.5*DELTAt*( bead[n+1].A*bead[n+1].m + bead[n].A*(bead[n].m - x0)) ) / ( 1.0 + 0.5*bead[n+1].A );
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

  //return -(x1-x0)*0.5*( (bead[n+1].A*(x1-bead[n+1].m)) + (bead[n].A*(x0-bead[n].m)) );
}

// =========================================
double pot(double x)
{
  // U = (x^2-1)^2
  //return 1.0 + x*x*(-2.0 + x*x);
  return (x*x-1.0)*(x*x-1.0);
}

// =========================================
//DeltaU: x m A
double DeltaU(double x, averages* bead, int n)
{
  return pot(x) - 0.5*bead[n].A*(x - bead[n].m)*(x - bead[n].m);
}

// =========================================
void writeConfig(averages* bead, int argc, char* argv[])
{
  int index;
  char filename[80];

  sprintf(filename, "%s-%s.dat", argv[3], argv[1]);
  FILE * pWriteStats;
  pWriteStats = fopen(filename, "w");

  for(index=0;index<NUMb;index++)
  {
    fprintf(pWriteStats, "%+.15e ", bead[index].m);
    fprintf(pWriteStats, "%+.15e ", bead[index].A);
    fprintf(pWriteStats, "%+.15e ", bead[index].xbar);
    fprintf(pWriteStats, "%+.15e ", bead[index].xxbar);
    fprintf(pWriteStats, "%+.15e ", bead[index].expVal[0]);
    fprintf(pWriteStats, "%+.15e ", bead[index].expVal[1]);
    fprintf(pWriteStats, "%+.15e ", bead[index].expVal[2]);
    fprintf(pWriteStats, "%+.15e ", bead[index].expVal[3]);
    fprintf(pWriteStats, "%+.15e ", bead[index].expVal[4]);
    fprintf(pWriteStats, "\n");
  }

  fclose(pWriteStats);
}
// =========================================
