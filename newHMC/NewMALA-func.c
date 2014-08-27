/*   NewMALA-func.c 
 *   Written Summer 2014 -- Patrick Malsom
 
function calls to the new implementation of the HMC algorithm

This is just a shared library that is linked with guided.py
*/

// ==================================================================================
// Definitions and Prototypes
// ==================================================================================
//STD Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Defines
#define NUML 99999
//
//Struct definitions
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

//Struct definitions
typedef struct _path
{
  double posOld;
  double posCur;
  double posNew;
} averages;

// ==================================================================================
// Function prototypes
// ==================================================================================
double Force(double x);
double ForcePrime(double x);
double Pot(double x);
double g(double xm1, double x0, double x1, parameters params, double invdt);

// ==================================================================================
double Pot(double x){
    return( 1.+ x*x*(-3.375+x * (1.6875 +x * (2.84766 +(-2.84766+0.711914 * x) * x))));
}

double Force(double x) {
  return (x * (6.75 + x * (-5.0625 + x * (-11.3906 + (14.2383 - 4.27148 * x) * x))));
}

double ForcePrime(double x){
    return (6.75 + x * (-10.125 + x * (-34.1719 + (56.9531 - 21.3574 * x) * x)));
}

double g(double xm1,double x0,double x1, parameters params, double invdt){
    double g1st=ForcePrime(x0)*Force((x0));
    double g2nd=-0.5*params.invdt*(Force((x1))-2.0*Force((x0))+Force((xm1)));
    double g3rd=-0.5*params.invdt*ForcePrime(x0)*(x1-2.0*x0+xm1);
    return g1st+g2nd+g3rd;
}
