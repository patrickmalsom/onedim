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

double gen_next_step(double x0,double random_gauss, parameters params){

  // save the initial point into two variables
  // these will be iterated upon in midpt and simpson
  double xsave = x0;
  double xnew = x0;

  switch(params.method) {
    case 0 :
      // return new step with leap frog quadrature
      xnew = x0 + params.dt * Force(x0) + params.noisePref * random_gauss;
      break;

    case 1 :
      // return new step with midpoint quadrature (implicit method)
      // guess the new step according to leap frog
      xnew = x0 + params.dt * Force(x0) + params.noisePref * random_gauss;

      // calculate new step iteratively with midpoint force
      while(abs(xnew-xsave) > 0.00000000001) {
        xsave=xnew;
        xnew = x0 + params.dt*midpt_Force(xnew,x0) + params.noisePref * random_gauss;
      }
      break;

    case 2 :
      // return new step with simpsons quadrature (implicit method)
      // guess the new step according to leap frog
      xnew = x0 + params.dt * Force(x0) + params.noisePref * random_gauss;

      //calculate new step iteratively with simsons force
      while(abs(xnew-xsave) > 0.00000000001) {
        xsave=xnew;
        xnew = x0 + params.dt*simpson_Force(xnew,x0) + params.noisePref * random_gauss;
      }
      break;
  }

  //return the final new step
  return(xnew);
}

double energy_drift(double x1, double x0, parameters params){

  double energy = 0.0;

  double Fx1 = Force(x1);
  double Fx0= Force(x0);

  switch(params.method) {
    case 0 :
      //leapfrog energy drift
      energy = 0.5*(x1-x0)*(Fx1+Fx0)+params.dt*0.25*(Fx1*Fx1-Fx0*Fx0)+Pot(x1)-Pot(x0);
      break;
    case 1 :
      //midpt energy drift
      energy = (x1-x0)*midpt_Force(x1,x0)+Pot(x1)-Pot(x0);
      break;
    case 2 :
      //simpson energy drift
      energy = (x1-x0)*simpson_Force(x1,x0)+Pot(x1)-Pot(x0);
      break;
  }

  return(energy);
}



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

  int loops = atoi(argv[7]);
  int suppress_print= atoi(argv[8]);

  params.noisePref = sqrt(2.0*params.eps*params.dt);


  //===============================================================
  // GNU Scientific Library Random Number Setup
  //===============================================================
  // Example shell command$ GSL_RNG_SEED=123 ./a.out
  const gsl_rng_type * RanNumType;
  gsl_rng *RanNumPointer; 
  gsl_rng_env_setup();
  RanNumType = gsl_rng_default;
  RanNumPointer= gsl_rng_alloc (RanNumType);

  //print parameters to stdout
  if(suppress_print == 1){
    printf("=======================================================\n");
    printf("dt:   %f\n",params.dt);
    printf("eps:  %f\n",params.eps);
    printf("x(0): %f\n",params.xstart);
    printf("num:  %i\n",params.num);
    printf("method {0:leapfrog,1:midpt,2:simpson}: %i\n",params.method);
    printf("MHMC {0:No,1:Yes}: %i\n",params.MHMC);
    printf("RNG: %s ", gsl_rng_name(RanNumPointer));
    printf("RNG Seed: %li \n", gsl_rng_default_seed);
    printf("=======================================================\n");
  }

  double inv_eps=1.0/params.eps;

  int acc;
  int Bs;

  double x1;
  double x0;

  int loop_iterator;
  int i;

  for(loop_iterator = 0; loop_iterator<loops; loop_iterator++){
    Bs = 0;
    acc = 0;
    x0 = params.xstart;

    for(i=0;i<params.num;i++){
      // generate the proposal move with params.method quadrature
      x1 = gen_next_step(x0,gsl_ran_gaussian(RanNumPointer,1),params);

      // perform the metropolis hastings test if params.MHMC==1
      if(params.MHMC == 1){ 
        //accept move if MHMC is satisfied
        if(exp( -energy_drift(x1,x0,params) * inv_eps) > gsl_rng_uniform(RanNumPointer)){
          x0 = x1;
          acc+=1;
        }
      }
      else{ 
        //always accept move  when MHMC is off (no rejections)
        x0=x1;
      }

      //calculate the broad well fraction (not averaged here)
      if(x0>0){
        Bs++;
      }
    }

    printf("%f\n",(float)Bs/params.num);
  }

  // free the GSL RNG memory pointers
  gsl_rng_free (RanNumPointer);

}
