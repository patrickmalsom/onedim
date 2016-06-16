#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// NOTE: The structs MUST EXACTLY MATCH the class in the python routine
// Stores many useful constants. initialized in the python code
typedef struct _parameters
{
  double dt;
  double eps;
  double noisePref;
  double xstart;
  int num;
  int method; //0:leapfrog 1:midpt 2:simpson
  int MHMC; //boolean for MHMC on:1 off:0
} parameters;

// Array struct for the trajectory
typedef struct _traj_array
{
  double grn;
  double pos;
  double rand;
} traj_array;

/* =====================================================================
 * Potential Definitions: Fatter-Skinnier Potential: (8 - 5*x)**8 * (2 + 5*x)**2 / 2**26
 * -----------------------------
 *   Pot:       U(x)          -> returns the potential at a position
 *   Force:     F(x)=-dU/dx   -> returns the force at a position
*/

// Calculate the fatter-skinny potential ( U(x) )
double Pot(double x){
  // Calculate: ((8 - 5 x)^8 (2. + 5 x)^2)/2^26
  // Horner form for the potential:
  return x*x*(x*(x*(x*(x*(x*(x*(x*(0.145519152283669*x - 1.74622982740402) + 8.96397978067398) - 25.331974029541) + 41.7232513427734) - 37.384033203125) + 10.68115234375) + 9.765625) - 7.8125) + 1.0; 


}

//Calculate the fatter skinny Force ( F = - dU(x)/dx )
double Force(double x){

  // Calculate: -250*x*(5*x-8)**7*(5*x+2)/2**26
  // Horner Form for the force
  return x*(x*(x*(x*(x*(x*(x*(x*(-1.45519152283669*x + 15.7160684466362) - 71.7118382453918) + 177.323818206787) - 250.339508056641) + 186.920166015625) - 42.724609375) - 29.296875) + 15.625);
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

double gen_leapfrog(double x0,double random_gauss, parameters params){
  // return new step 
  return x0 + params.dt * Force(x0) + params.noisePref * random_gauss;
}

double gen_midpt(double x0,double random_gauss, parameters params){
  // save the initial point
  double xsave = x0;
  // guess the new step according to leap frog
  double xguess = x0 + params.dt * Force(x0) + params.noisePref * random_gauss;

  //calculate new step iteratively with midpoint force
  while(abs(xguess-xsave) > 0.00000000001) {
    xsave=xguess;
    xguess = x0 + params.dt*midpt_Force(xguess,x0) + params.noisePref * random_gauss;
  }
  return xguess;
}

double gen_simpson(double x0,double random_gauss, parameters params){
  // save the initial point
  double xsave = x0;
  // guess the new step according to leap frog
  double xguess = x0 + params.dt * Force(x0) + params.noisePref * random_gauss;

  //calculate new step iteratively with simsons force
  while(abs(xguess-xsave) > 0.00000000001) {
    xsave=xguess;
    xguess = x0 + params.dt*simpson_Force(xguess,x0) + params.noisePref * random_gauss;
  }
  return xguess;
}

int create_trajectory(traj_array* traj, parameters params){
  // Create the whole trajectory and keep track of B(s)
  int i;
  int MHMC_acc = 0;

  double inv_eps=1.0/params.eps;

  double MHMC_criteria;
  double proposed_move;

  //save value of the broad well fraction
  traj[0].pos=params.xstart;

  // leapfrog trajectory without MHMC
  if(params.method==0 && params.MHMC==0){
    for(i=0;i<params.num-1;i++){
      traj[i+1].pos = gen_leapfrog(traj[i].pos,traj[i].grn,params);
    }
  }

  //TODO THIS IS NOT CORRECT!
  // leapfrog trajectory with MHMC
  if(params.method==0 && params.MHMC==1){
    // calculate 1/(2*eps) for use later


    for(i=0;i<params.num-1;i++){
      // generate the proposal move
      proposed_move = gen_leapfrog(traj[i].pos,traj[i].grn,params);

      // perform the metropolis hastings test
      MHMC_criteria = exp( -( Pot(proposed_move) - Pot(traj[i].pos) )*inv_eps);
      if(MHMC_criteria > traj[i].rand){
        traj[i+1].pos=proposed_move;
        MHMC_acc+=1;
      }
      else{
        traj[i+1].pos=traj[i].pos;
      }
    }
  }

  // midpt trajectory without MHMC
  if(params.method==1 && params.MHMC==0){
    for(i=0;i<params.num-1;i++){
      traj[i+1].pos = gen_midpt(traj[i].pos,traj[i].grn,params);
    }
  }

  // leapfrog trajectory without MHMC
  if(params.method==2 && params.MHMC==0){
    for(i=0;i<params.num-1;i++){
      traj[i+1].pos = gen_simpson(traj[i].pos,traj[i].grn,params);
    }
  }
  return(MHMC_acc);
}

