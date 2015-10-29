
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
 * Potential Definitions: Fatter-Skinnier Potential: ((8 - 5 x)^8 (2. + 5 x)^2)/2**26
 * -----------------------------
 *   Pot:       U(x)          -> returns the potential at a position
 *   Force:     F(x)=-dU/dx   -> returns the force at a position
 *   ForcePr:   F'(x)=dF/dx   -> returns the first deriv of force
 *   ForcePrPr: F''(x)=dF'/dx -> returns the second deriv of force
*/

double Pot(double x){
  return 1.0000000000000002 + x*x*(-7.812499999999998 + x*(9.765625 + x*(10.68115234375 + x*(-37.384033203125 + x*(41.72325134277344 + x*(-25.331974029541016 + x*(8.96397978067398 + (-1.7462298274040222 + 0.14551915228366852*x)*x)))))));
}

double Force(double x){
  return x*(15.625 + x*(-29.296875 + x*(-42.724609375 + x*(186.920166015625 + x*(-250.33950805664065 + x*(177.3238182067871 + x*(-71.71183824539185 + (15.716068446636202 - 1.4551915228366852*x)*x)))))));
}

double ForcePrime(double x){
  return 15.625 + x*(-58.59375 + x*(-128.173828125 + x*(747.6806640625 + x*(-1251.6975402832031 + x*(1063.9429092407227 + x*(-501.9828677177429 + (125.7285475730896 - 13.096723705530167*x)*x))))));
}

double ForceDoublePrime(double x){
  return -58.59375 + x*(-256.34765625 + x*(2243.0419921875 + x*(-5006.7901611328125 + x*(5319.714546203613 + x*(-3011.8972063064575 + (880.0998330116272 - 104.77378964424133*x)*x)))));
}

