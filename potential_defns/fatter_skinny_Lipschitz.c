
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
 * Potential Definitions: Fatter-Skinnier Potential: 
 * Globally Lipshitz: 2^-26 ((8 - 5 x)^8 (2. + 5 x)^2)/(1 + (x/4)^8)
 * -----------------------------
 *   Pot:       U(x)          -> returns the potential at a position
 *   Force:     F(x)=-dU/dx   -> returns the force at a position
 *   ForcePr:   F'(x)=dF/dx   -> returns the first deriv of force
 *   ForcePrPr: F''(x)=dF'/dx -> returns the second deriv of force
*/

double Pot(double x){
  double x2=x*x;
  double x8=x2*x2*x2*x2;
  return (4.398046511104e12 + x2*(-3.4359738368e13 + x*(4.294967296e13 + x*(4.697620479999999e13 + x*(-1.644167168e14 + x*(1.8350080000000003e14 + x*(-1.1141120000000002e14 + x*(3.9424e13 + x*(-7.68e12 + 6.4e11*x)))))))))/(67108864*(65536. + x8));
}

double Force(double x){
  double x2=x*x;
  double x8=x2*x2*x2*x2;
  return (x*(6.7108864e10 + x*(-1.2582912e11 + x*(-1.835008e11 + x*(8.02816e11 + x*(-1.0752e12 + x*(7.616e11 + x*(-3.07999475712e11 + x*(6.75e10 + x*(-6.253072e9 + x*(3.2e6 + x*(2.8e6 + x*(-7.35e6 + x*(5.46875e6 + x*(-1.66015625e6 + (114440.91796875 - 19073.486328125*x)*x2)))))))))))))))/(4.294967296e9 + x8*(131072. + x8));
}

double ForcePrime(double x){
  double x2=x*x;
  double x8=x2*x2*x2*x2;
  return (4.398046511104e15 + x*(-1.649267441664e16 + x*(-3.60777252864e16 + x*(2.10453397504e17 + x*(-3.52321536e17 + x*(2.994733056e17 + x*(-1.4129537548183142e17 + x*(3.538944e16 + x*(-3.6892185722880005e15 + x*(3.85875968e12 + x*(4.4040192e12 + x*(-1.54140672e13 + x*(1.6486399999999998e13 + x*(-9.1392e12 + x*(2.771995281408e12 + x*(-4.2e11 + x*(2.2521503999999996e10 + x*(-1.92e7 + x*(-1.4000000000000002e7 + x*(2.94e7 + x*(-1.640625e7 + x*(3.320312499999985e6 - 19073.486328125*x*x2))))))))))))))))))))))/(2.81474976710656e14 + x8*(1.2884901888e10 + x8*(196608. + x8)));
}

double ForceDoublePrime(double x){
  double x2=x*x;
  double x8=x2*x2*x2*x2;
  return (-1.0808639105689189e21 + x*(-4.728779608739019e21 + x*(4.137682157646642e22 + x*(-9.235897673318398e22 + x*(9.813141277900798e22 + x*(-5.5559602365463825e22 + x*(1.623497637888e22 + x*(-1.934318579943997e21 + x*(2.65532058107904e18 + x*(3.6799279792127995e18 + x*(-1.55314607357952e19 + x*(2.00118632448e19 + x*(-1.3476298751999998e19 + x*(5.086633517345931e18 + x*(-1.0144972799999999e18 + x*(8.2643005734912e16 + x*(-7.927234559999998e13 + x*(-7.817134079999998e13 + x*(2.369912832e14 + x*(-2.193408e14 + x*(1.0510079999999998e14 + x*(-2.771995281408e13 + x*(3.78e12 + x*(-2.1017203199999954e11 + x*(1.344e8 + x*(8.399999999999999e7 + x*(-1.4699999999999997e8 + (6.5624999999999985e7 - 9.9609375e6*x)*x)))))))))))))))))))))))))))/(1.8446744073709552e19 + x8*(1.125899906842624e15 + x8*(2.5769803776e10 + x8*(262144. + x8))));
}

