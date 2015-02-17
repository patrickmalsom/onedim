
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
 * Potential Definitions: Fat-Skinny Potential: (3x-4)^4*(3x+2)^2/1024
 * -----------------------------
 *   Force:     F(x)=-dU/dx   -> returns the force at a position
 *   ForcePr:   F'(x)=dF/dx   -> returns the first deriv of force
 *   ForcePrPr: F''(x)=dF'/dx -> returns the second deriv of force
*/
double Force(double x){
  return x*(6.75 + x*(-5.0625 + x*(-11.390625 + (14.23828125 - 4.271484375*x)*x)));
}
double ForcePrime(double x){
  return 6.75 + x*(-10.125 + x*(-34.171875 + (56.953125 - 21.357421875*x)*x));
}
double ForceDoublePrime(double x){
  return -10.125 + x*(-68.34375 + (170.859375 - 85.4296875*x)*x);
}
