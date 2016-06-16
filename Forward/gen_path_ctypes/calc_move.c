#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//GNU Scientific Libraries
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

typedef struct _parameters{
  double dt;
  double eps;
  double noisePref;
  double xstart;
  int num;
  int method; //0:leapfrog 1:midpt 2:simpson
  int MHMC; //boolean for MHMC on:1 off:0
} parameters;


double Pot(double x){
  // Calculate the fatter-skinny potential ( U(x) )
  // Horner form of: ((8 - 5 x)^8 (2. + 5 x)^2)/2^26
  return x*x*(x*(x*(x*(x*(x*(x*(x*(0.145519152283669*x - 1.74622982740402) + 8.96397978067398) - 25.331974029541) + 41.7232513427734) - 37.384033203125) + 10.68115234375) + 9.765625) - 7.8125) + 1.0; 
}


double Force(double x){
  // Calculate the fatter skinny Force ( F = - dU(x)/dx )
  // Horner form of: -250*x*(5*x-8)**7*(5*x+2)/2**26
  return x*(x*(x*(x*(x*(x*(x*(x*(-1.45519152283669*x + 15.7160684466362) - 71.7118382453918) + 177.323818206787) - 250.339508056641) + 186.920166015625) - 42.724609375) - 29.296875) + 15.625);
}


double midpt_Force(double x1,double x0){
  // return the force at the mid point
  // F_w(x1,x0) = F( (x1+x0)/2 )
  return Force(0.5*(x1+x0));
}


double simpson_Force(double x1,double x0){
  // return the force according to simpsons method
  // F_w(x1,x0) = ( F(x1) + 4*F((x1+x0)/2) + F(x0) )/6
  return (Force(x1) + 4.0*Force(0.5*(x1+x0)) + Force(x0))*0.16666666666666666666;
}


double gen_leapfrog(double x0,double random_gauss, parameters params){
  // return new step with leap frog quadrature
  return x0 + params.dt * Force(x0) + params.noisePref * random_gauss;
}

double gen_midpt(double x0,double random_gauss, parameters params){
  // return new step with midpoint quadrature (implicit method)
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
  // return new step with simpsons quadrature (implicit method)
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

/*
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
*/









int main(int argc, char *argv[])
{

  //Declare params
  parameters params;
  //Initialize the params struct
  // argument order: dt eps noisePref xstart num method MHMC
  // the order needs to be changed in the bash script if it is changed here!
  params.dt = atof(argv[1]);
  params.eps = atof(argv[2]);
  params.xstart = atof(argv[3]);
  params.num = atoi(argv[4]);
  params.method = atoi(argv[5]);
  params.MHMC = atoi(argv[6]);

  params.noisePref = sqrt(2.0*params.eps*params.dt);

  //===============================================================
  // GNU Scientific Library Random Number Setup
  //===============================================================
  // Example shell command$ GSL_RNG_SEED=123 ./a.out
  printf("=======================================================\n");
  const gsl_rng_type * RanNumType;
  gsl_rng *RanNumPointer; 
  gsl_rng_env_setup();
  RanNumType = gsl_rng_default;
  RanNumPointer= gsl_rng_alloc (RanNumType);
  printf("Random Number Generator Type: %s \n", gsl_rng_name(RanNumPointer));
  printf("RNG Seed: %li \n", gsl_rng_default_seed);
  printf("=======================================================\n");

  //printf("%f\n",gsl_rng_uniform(RanNumPointer));
  //printf("%f\n",gsl_ran_gaussian(RanNumPointer,1));

  double x1;
  double x0 = params.xstart;
  int acc = 0;
  int Bs = 0;

  int i;
  double inv_eps=1.0/params.eps;
  double Fx1;
  double Fx0;
  double energy;

  for(i=0;i<params.num;i++){
    // generate the proposal move
    x1 = gen_leapfrog(x0,gsl_ran_gaussian(RanNumPointer,1),params);

    // perform the metropolis hastings test
    Fx1= Force(x1);
    Fx0= Force(x0);
    energy = 0.5*(x1-x0)*(Fx1+Fx0)+params.dt*0.25*(Fx1*Fx1-Fx0*Fx0)+Pot(x1)-Pot(x0);

    if(exp( -energy * inv_eps) > gsl_rng_uniform(RanNumPointer)){
      x0 = x1;
      acc+=1;
    }
    if(x0>0){
      Bs++;
    }
  }

  printf("%f\n",(float)Bs/params.num);
  gsl_rng_free (RanNumPointer);

}
