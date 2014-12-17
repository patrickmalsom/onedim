/*   ItoHMC.c 
     Written Winter 2014 -- Patrick Malsom
 
C library for the Ito HMC algorithm
Functions used to generate the SPDE and the MD steps

This is a shared library that is linked to python (ItoHMC.py)
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
} parameters;

// Path Struct
// stores an array of useful quantities for the SPDE and MD simulation
typedef struct _path
{
  double pos;
  double randlist;
  double bb;
  double F;
  double Fp;
  double Fpp;
  double G;
  double gradG;
  double LinvG;
} averages;

// ==================================================================================
// Functions
// ==================================================================================

// ==================================================================================
//fat skinny
double Pot(double x){
  // U(x): returns the potential at a position
  return 1. + x*x*(-3.375 + x*(1.6875 + x*(2.84765625 + (-2.84765625 + 0.7119140625*x)*x)));
}
double Force(double x) {
  // F(x)=-dU/dx: returns the force at a position
  return x*(6.75 + x*(-5.0625 + x*(-11.390625 + (14.23828125 - 4.271484375*x)*x)));
}
double ForcePrime(double x){
  // F'(x)=dF/dx: returns the force at a position
  return 6.75 + x*(-10.125 + x*(-34.171875 + (56.953125 - 21.357421875*x)*x));
}
double ForceDoublePrime(double x){
  // F'(x)=dF/dx: returns the force at a position
  return -10.125 + x*(-68.34375 + (170.859375 - 85.4296875*x)*x);
}

//============================================
void calcPotentials(averages* path, parameters params)
{
  int i;

  for(i=0;i<params.NumB;i++){

    path[i].F=Force(path[i].pos);
    path[i].Fp=ForcePrime(path[i].pos);
    path[i].Fpp=ForceDoublePrime(path[i].pos);

    path[i].G=0.5*path[i].F*path[i].F + params.eps*path[i].Fp;
    path[i].gradG=path[i].F*path[i].Fp + params.eps*path[i].Fpp;
   
  }
}

//============================================
void LInverse(averages* path, parameters params)
{
  double lasti;
  int n;
  int Numl=params.NumB-2;
  double vecdg[Numl];
  double veci0[Numl];
  double veci1[Numl];

  for(n=0;n<params.NumB-2;n++){
    vecdg[n]=path[n+1].gradG*params.deltat*params.deltat;
  }
  veci0[0]=vecdg[0];
  veci1[0]=vecdg[0];
  //this for loop is recursive!!!!
  for(n=1;n<Numl;n++){
    veci0[n]=veci0[n-1]+vecdg[n];
    veci1[n]=veci1[n-1]+((double)(n+1))*vecdg[n];
  }

  lasti=veci0[Numl-1]-(path[params.NumB-1].pos-path[0].pos+veci1[Numl-1])/((double)(Numl+1));

  for(n=0;n<Numl;n++){
    path[n+1].LinvG=path[0].pos+((double)(n+1))*(veci0[n]-lasti) - veci1[n];
  }
  path[0].LinvG=path[0].pos;
  path[params.NumB-1].LinvG=path[params.NumB-1].pos;

}

// ==================================================================================
void calcSPDEpos(averages* path0, averages* path1, parameters params){
  // calculate the right hand side vector for the SPDE
  // this is x^{(1)} vector in the notes

  int i;
  double h=sqrt(2.0*params.deltat);
  double si=(4.0*h)/(4.0+h*h);
  double co=(4.0-h*h)/(4.0+h*h);
  double noisePref=sqrt(2.0*params.eps*params.deltat);

  for(i=1;i<params.NumB-1;i++){
  path1[i].pos=si*h*0.5*path0[i].LinvG + si*noisePref*path0[i].randlist + co * path0[i].pos;
  }

  //set boundary conditions
  path1[0].pos=path0[0].pos;
  path1[params.NumB-1].pos=path0[params.NumB-1].pos;

}

// ==================================================================================
void calcMDpos(averages* path0, averages* path1, averages* path2, parameters params){
  // calculate the right hand side vector for the MD steps
  // this is x^{(2)} vector in the notes

  int i;
  double h=sqrt(2.0*params.deltat);
  double si=(4.0*h)/(4.0+h*h);
  double co=(4.0-h*h)/(4.0+h*h);

  for(i=1;i<params.NumB-1;i++){
  path2[i].pos=(path1[i].pos-path0[i].pos) + (2.0*co - 1.0)*path1[i].pos + si*h*path1[i].LinvG;
  }

  //set boundary conditions
  path2[0].pos=path1[0].pos;
  path2[params.NumB-1].pos=path1[params.NumB-1].pos;

}

// ==================================================================================
double calcEnergyChange(averages* path0, averages* path1, parameters params){

  int i;

  double h=sqrt(2.0*params.deltat);
  double cot=(4.0l-h*h)/(4.0l*h);
  double csc=(4.0l+h*h)/(4.0l*h);

  double lambda1=0.0;
  double lambda2=0.0;
  // tempsum is a temporary double for storage of the integrals in the lambdas
  // it should be set to zero at the beginning of each integration evaluation.
  double tempsum;

  // first integral in Lambda1
  tempsum=0.0;
  for(i=0;i<params.NumB;i++){
    tempsum+= (path1[i].pos *path0[i].gradG - path0[i].pos*path1[i].gradG);
  }
  lambda1-= h/(4.0*params.eps)*csc * params.deltat*tempsum;

  // second integral in Lambda1
  tempsum=0.0;
  for(i=0;i<params.NumB;i++){
    tempsum+= (path0[i].pos*path0[i].gradG - path1[i].pos*path1[i].gradG);
  }
  lambda1+=h/(4.0*params.eps)*cot * params.deltat*tempsum;

  // third integral in Lambda1
  tempsum=0.0;
  for(i=0;i<params.NumB;i++){
    tempsum+= ( path0[i].gradG*path0[i].LinvG - path1[i].gradG*path1[i].LinvG );
  }
  lambda1+=h*h/(16.0*params.eps) * params.deltat*tempsum;

  // integral in Lambda2
  tempsum=0.0;
  for(i=0;i<params.NumB;i++){
    tempsum+= (path1[i].G - path0[i].G);
  }
  lambda2+=0.5/params.eps * params.deltat*tempsum;

  //return the change in energy for a path step
  return( lambda1+lambda2 );

}


//============================================
void generateBB(averages* path, parameters params){
//generate a Brownian bridge and store to BrownianBridge
  int n;
  double xn;
  double doubleNUMu=((double)(params.NumB-1));

  double endPtCorr;
  double sum, term, term0;
  double alpha;

  double sqdu=sqrt(2.0*params.eps*params.deltat);
  path[0].bb=0.0l;
  for(n=1;n<params.NumB;n++)
  {
    path[n].bb=path[n-1].bb+sqdu*path[n-1].randlist;
  }

  xn=path[params.NumB-1].bb/((double)(params.NumB-1));
  for(n=1;n<params.NumB-1;n++)
  {
    path[n].bb-=((double)(n))*xn;
  }
  path[params.NumB-1].bb=0.0l;

  //normalize the brownian bridge to have accurate quad var
  endPtCorr= (path[0].bb-path[params.NumB-1].bb)*(path[0].bb-path[params.NumB-1].bb)/doubleNUMu;
  sum=0.0l;
  for(n=0;n<params.NumB-1;n++){
    sum+= (path[n].bb-path[n+1].bb)*(path[n].bb-path[n+1].bb);
  }
  alpha=sqrt((doubleNUMu*params.deltat*2.0*params.eps-endPtCorr)/(sum-endPtCorr));
  term=(1.0l-alpha)*(path[params.NumB-1].bb-path[0].bb)/doubleNUMu;
  term0=(1.0l-alpha)*path[0].bb;
  //term and term0 have a subtraction of roughly equal numbers and thus is not very accurate
  // alpha is ~1 with an error of 10^-4 or 5 for sample configs. This makes the routine
  //nondeterministic between Fortran and C
  for(n=1;n<params.NumB-1;n++){
    path[n].bb=alpha*path[n].bb+term0+((double)(n-1))*term;
  }
}

// ==================================================================================
double quadVar(averages* path, parameters params){
  // fill out the Hessian variables for the structure
  // initial params must have the positions filled
  double qv=0.0;
  int i;

  for(i=0;i<params.NumB-1;i++){
    qv+=(path[i+1].pos-path[i].pos)*(path[i+1].pos-path[i].pos);
  }
  return(qv/(2.0*params.eps*params.deltat*(params.NumB-1)));
}

