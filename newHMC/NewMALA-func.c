/*   NewMALA-func.c 
     Written Summer 2014 -- Patrick Malsom
 
C library for the finite time double HMC algorithm
Functions used to generate the SPDE and the MD steps

This is a shared library that is linked to python (NewMALA.py)
*/

//STD Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ==================================================================================
// Data Structures
// ==================================================================================
// IMPORTANT: These structs must EXACTLY MATCH the structs defined in the pyton code!

// Parameters Struct 
// Stores many useful constants that are defined in the python code
// ( see python code for comments)
typedef struct _parameters
{
  double deltat;
  double invdt;
  double eps;
  double deltatau;
  double noisePref;
  double r;
  int NumB;
  int NumU;
  int NumL;
  double xPlus;
  double xMinus;
} parameters;

// Path Struct
// stores an array of useful quantities for the SPDE and MD simulation
//     pos: array of positions 
//     randlist: random gaussian numbers(0,1). should be passed from python 
//               (TODO: read about the numpy RNG)
//     Force/ForcePrime: array to store the forces along the path
//     dg: dG/dx array
//     rhs: RHS of the SPDE/MD eqns (different eqns for SPDE and MD)
typedef struct _path
{
  double pos;
  double randlist;
  double Force;
  double ForcePrime;
  double dg;
  double rhs;
} averages;

// ==================================================================================
// Function prototypes
// ==================================================================================

double Pot(double x);
// U(x): returns the potential at a position

double Force(double x);
// F(x)=-dU/dx: returns the force at a position

double ForcePrime(double x);
// F'(x)=dF/dx: returns the force at a position

void calcForces(averages* path, parameters params);
void calcForcePrime(averages* path, parameters params);
void calcdg(averages* path, parameters params);
void calcSPDErhs(averages* path, parameters params);
void calcStateSPDE(averages* path, parameters params);

void GaussElim(averages* path0, averages* path1, parameters params);
// Gaussian elimination for computing L.x=b where L is the second deriv matrix
//   L matrix: mainDiag: (1+2r) ; upper(lower)Diag: (-r)
//   input: the filled path0 state (
//   output: new path1.pos  


// ==================================================================================
// Functions
// ==================================================================================

double Pot(double x){
  // U(x): returns the potential at a position
  return( 1.+ x*x*(-3.375+x * (1.6875 +x * (2.84766 +(-2.84766+0.711914 * x) * x))));
}

double Force(double x) {
  // F(x)=-dU/dx: returns the force at a position
  return (x * (6.75 + x * (-5.0625 + x * (-11.3906 + (14.2383 - 4.27148 * x) * x))));
}

double ForcePrime(double x){
  // F'(x)=dF/dx: returns the force at a position
  return (6.75 + x * (-10.125 + x * (-34.1719 + (56.9531 - 21.3574 * x) * x)));
}

// ==================================================================================
void calcForces(averages* path, parameters params){
// fill out the Force variables for the struce
// initial params must have the positions filled
  double x;
  int i;

  for(i=0;i<params.NumB;i++){
    x = path[i].pos;
    path[i].Force= (x * (6.75 + x * (-5.0625 + x * (-11.3906 + (14.2383 - 4.27148 * x) * x))));
  }
}

// ==================================================================================
void calcForcePrimes(averages* path, parameters params){
// fill out the ForcePrime variables for the struce
// initial params must have the positions filled
  double x;
  int i;

  for(i=0;i<params.NumB;i++){
    x = path[i].pos;
    path[i].ForcePrime= (6.75 + x * (-10.125 + x * (-34.1719 + (56.9531 - 21.3574 * x) * x)));
  }
}


// ==================================================================================
void calcdg(averages* path, parameters params){
// calculate dG/dx
    int i;
    double g1st, g2nd, g3rd;
    for(i=1;i<params.NumB-1;i++){
      g1st=path[i].ForcePrime*path[i].Force;
      g2nd=-0.5*params.invdt*(path[i+1].Force-2.0*path[i].Force+path[i-1].Force);
      g3rd=-0.5*params.invdt*path[i].ForcePrime*(path[i+1].pos-2.0*path[i].pos+path[i-1].pos);
      path[i].dg=params.deltatau*(g1st+g2nd+g3rd);
    }
    //printf("dg=%+0.15e\n",params.noisePref*path[1].randlist);
}

// ==================================================================================
void calcSPDErhs(averages* path, parameters params){
  // x vector
  // note that the i_th position is (i+1)
  int i;
  double temp1, temp2, temp3;

  for(i=1;i<params.NumB-1;i++){
    temp1 = params.r*path[i-1].pos + (1.-2.*params.r)*path[i].pos + params.r*path[i+1].pos;
    temp2 = path[i].dg;
    temp3 = params.noisePref*path[i].randlist;
    path[i].rhs=temp1-temp2+temp3;
  }
  path[1].rhs+=params.r*path[0].pos;
  path[params.NumB-2].rhs+=params.r*path[params.NumB-1].pos;

}

// ==================================================================================
void calcStateSPDE(averages* path, parameters params){
// Function to calculate all of the structure parameters 
// using the pos variables and the randlist

  calcForces(path, params);
  calcForcePrimes(path, params);
  calcdg(path, params);
  calcSPDErhs(path, params);
}

/*
// ==================================================================================
void calcStateMD(averages* path, parameters params){
// Function to calculate all of the structure parameters 
// using the pos variables

  calcForces(path,params);
  calcForcePrimes(path,params);
  calcdg(path,params);
  calcrhsMD(path,params);
}

// ==================================================================================
void calcMDrhs(averages* path0, averages* path1, averages* path2, parameters params){
}

// ==================================================================================
double calcEnergyChange(averages* path0, averages* path1, parameters params){
}
*/


// ==================================================================================
// Gaussian elimination for computing L.x=b where
//   L is tridiagonal matrix
//     mainDiag: 1+2r
//     upper(lower)Diag: -r
//
//   input: the filled path0 state 
//   output: new path1.pos  
void GaussElim(averages* path0, averages* path1, parameters params){

  int NumL = params.NumL;
  int i;
  double temp;
  
  // declare the L matrix arrays
  double al[NumL-1];
  double am[NumL];
  double au[NumL-1];

  //initialize the tridiagonal arrays
  for(i=0;i<NumL-1;i++){ 
    al[i]=-params.r;
    am[i]=(1.0+2.0*params.r);
    au[i]=-params.r;
  }
  am[NumL-1]=(1.0+2.0*params.r);


  // Gaussian Elimination
  for(i=0;i<NumL-1;i++){
    temp=-al[i]/am[i];
    am[i+1]+=au[i]*temp;
    path0[i+2].rhs+=path0[i+1].rhs*temp;
  }

  // Back substitution
  for(i=0;i<NumL-1;i++){
    temp=-au[NumL-i-2]/am[NumL-i-1];
    path0[NumL-i-1].rhs+=path0[NumL-i].rhs*temp;
  }
  
  // Divide by main diagonal
  for(i=0;i<NumL;i++){
    path0[i+1].rhs=path0[i+1].rhs/am[i];
  }

  // Now the rhs vector transformed to the new position vector
  // save the new positions to the path1 positions
  for(i=1;i<params.NumB-1;i++){
    path1[i].pos=path0[i].rhs;

  // set the boundary conditions on the new path positions
  // BCs are not included in the Gaussan elimination
  path1[0].pos=path0[0].pos;
  path1[params.NumB-1].pos=path0[params.NumB-1].pos;
  }
}

// ==================================================================================

