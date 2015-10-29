
//STD Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// include the struct for HMC simulation
#include "../onedimHMC_struct.h"
#include "../onedimHMC_forces.h"

//Include the HMC library
#include "../onedimHMC.c"


/* =====================================================================
 * Potential Definitions: quartic Potential: x^8
 * -----------------------------
 *   Force:     F(x)=-dU/dx   -> returns the force at a position
 *   ForcePr:   F'(x)=dF/dx   -> returns the first deriv of force
 *   ForcePrPr: F''(x)=dF'/dx -> returns the second deriv of force
*/
double Pot(double x){
  double xx = x*x;
  return xx*xx*xx*xx;
}
double Force(double x){
  double xx = x*x;
  return -8.0*xx*xx*xx*x;
}
double ForcePrime(double x){
  double xx = x*x;
  return -56.0*xx*xx*xx;
}
double ForceDoublePrime(double x){
  double xx = x*x;
  return -336.0*xx*xx*x;
}
