
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
  return x*x*(x*(x*(x*(x*(x*(x*(x*(0.145519152283669*x - 1.74622982740402) + 8.96397978067398) - 25.331974029541) + 41.7232513427734) - 37.384033203125) + 10.68115234375) + 9.765625) - 7.8125) + 1.0;
}

double Force(double x){
  return x*(x*(x*(x*(x*(x*(x*(x*(-1.45519152283669*x + 15.7160684466362) - 71.7118382453918) + 177.323818206787) - 250.339508056641) + 186.920166015625) - 42.724609375) - 29.296875) + 15.625);
}

double ForcePrime(double x){
  return x*(x*(x*(x*(x*(x*(x*(-13.0967237055302*x + 125.72854757309) - 501.982867717743) + 1063.94290924072) - 1251.6975402832) + 747.6806640625) - 128.173828125) - 58.59375) + 15.625;
}

double ForceDoublePrime(double x){
  return x*(x*(x*(x*(x*(x*(-104.773789644241*x + 880.099833011627) - 3011.89720630646) + 5319.71454620361) - 5006.79016113281) + 2243.0419921875) - 256.34765625) - 58.59375;
}
