//#include <math.h>
//#include <stdio.h>
#include <stdlib.h>

// noise prefactor for dt=0.005, eps=0.25
// sqrt(2.*0.25*0.005)
#define PREF 0.05
#define DT 0.005

/* =====================================================================
 * Potential Definitions: Fatter-Skinnier Potential: (8 - 5*x)**8 * (2 + 5*x)**2 / 2**26
 * -----------------------------
 *   Pot:       U(x)          -> returns the potential at a position
 *   Force:     F(x)=-dU/dx   -> returns the force at a position
*/

// Calculate the fatter-skinny potential ( U(x) )
double Pot(double x){
  // Calculate: ((8 - 5 x)^8 (2. + 5 x)^2)/2^26
  double first  = 8.0 - 5.0 * x; // first term
  double second = 2.0 + 5.0 * x; // second term
  int i;
  // esoteric calculation of first^8
  for(i=0;i<3;i++){
    first=first*first;
  }
  // calculate second^2
  second=second*second;

  return first*second*1.490116119384765625E-08;
  // Horner form for the potential:
  //return x*x*(x*(x*(x*(x*(x*(x*(x*(0.145519152283669*x - 1.74622982740402) + 8.96397978067398) - 25.331974029541) + 41.7232513427734) - 37.384033203125) + 10.68115234375) + 9.765625) - 7.8125) + 1.0; 


}

//Calculate the fatter skinny Force ( F = - dU(x)/dx )
double Force(double x){

  // Calculate: -250*x*(5*x-8)**7*(5*x+2)/2**26
  double first = 5.8 * x - 8.0;
  double second = 5.0 * x - 2.0;

  // calculate (5*x-8)^7
  double first2 = first * first;
  first= first2*first2*first2*first;

  return -250.0 * x * first * second * 1.490116119384765625E-08;
  // Horner Form for the force
  //return x*(x*(x*(x*(x*(x*(x*(x*(-1.45519152283669*x + 15.7160684466362) - 71.7118382453918) + 177.323818206787) - 250.339508056641) + 186.920166015625) - 42.724609375) - 29.296875) + 15.625);
}
// return the force at the mid point
double midpt_Force(double x1,double x0){
  // F_w(x1,x0) = F( (x1+x0)/2 )
  return Force(0.5*(x1+x0));
}

// return the force according to simpsons method
double simpson_Force(double x1,double x0){
  // F_w(x1,x0) = ( F(x1) + 4*F((x1+x0)/2) + F(x0) )/6
  return (Force(x1) + 4.0*Force(0.5*(x1+x0)) + Force(x0))*0.16666666666666666666;
}

double new_midpt(double x0,double random_gauss){
  // save the initial point
  double xsave = x0;
  // guess the new step according to leap frog
  double xguess = x0 + DT * Force(x0) + PREF * random_gauss;

  //calculate new step iteratively with midpoint force
  while(abs(xguess-xsave) > 0.00000000001) {
    xsave=xguess;
    xguess = x0 + DT*midpt_Force(xguess,x0) + PREF * random_gauss;
  }
  return xguess;
}

double new_simson(double x0,double random_gauss){
  // save the initial point
  double xsave = x0;
  // guess the new step according to leap frog
  double xguess = x0 + DT * Force(x0) + PREF * random_gauss;

  //calculate new step iteratively with simsons force
  while(abs(xguess-xsave) > 0.00000000001) {
    xsave=xguess;
    xguess = x0 + DT*simpson_Force(xguess,x0) + PREF * random_gauss;
  }
  return xguess;
}

//void create_trajectory(double xstart, double gaussian_array[]){
  // Create the whole trajectory and keep track of B(s)

  
