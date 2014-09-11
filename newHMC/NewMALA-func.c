/*   NewMALA-func.c 
     Written Summer 2014 -- Patrick Malsom
 
C library for the finite time double HMC algorithm
Functions used to generate the SPDE and the MD steps

This is a shared library that is linked to python (NewMALA.py)
*/

//STD Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ==================================================================================
// Data Structures
// ==================================================================================
// IMPORTANT: These structs must EXACTLY MATCH the structs defined in the pyton code!

// Parameters Struct 
// Stores many useful constants that are defined in the python code
// ( see python code for comments)
typedef struct _parameters
{
  double deltat;
  double invdt;
  double eps;
  double deltatau;
  double noisePref;
  double r;
  int NumB;
} parameters;

// Path Struct
// stores an array of useful quantities for the SPDE and MD simulation
//     pos: array of positions 
//     randlist: random gaussian numbers(0,1). should be passed from python 
//               (TODO: read about the numpy RNG)
//     Force/Hessian: array to store the forces/Hessian along the path
//     dg: dG/dx array
//     rhs: RHS of the SPDE/MD eqns (different eqns for SPDE and MD)
typedef struct _path
{
  double pos;
  double randlist;
  double Force;
  double Hessian;
  double deltae;
  double dg;
  double Phi;
  double rhs;
} averages;

// ==================================================================================
// Functions
// ==================================================================================

double Pot(double x){
  // U(x): returns the potential at a position
  return( 1.+ x*x*(-3.375+x * (1.6875 +x * (2.84766 +(-2.84766+0.711914 * x) * x))));
}

double Force(double x) {
  // F(x)=-dU/dx: returns the force at a position
  return (x * (6.75 + x * (-5.0625 + x * (-11.3906 + (14.2383 - 4.27148 * x) * x))));
}

double Hessian(double x){
  // F'(x)=dF/dx: returns the force at a position
  return -(6.75 + x * (-10.125 + x * (-34.1719 + (56.9531 - 21.3574 * x) * x)));
}


// ==================================================================================
void calcForces(averages* path, parameters params){
  // fill out the Force variables for the struce
  // initial params must have the positions filled
  double x;
  int i;

  for(i=0;i<params.NumB;i++){
    x = path[i].pos;
    //path[i].Force= (x * (6.75 + x * (-5.0625 + x * (-11.3906 + (14.2383 - 4.27148 * x) * x))));
    path[i].Force= x*(6.75 + x*(-5.0625 + x*(-11.390625 + (14.23828125 - 4.271484375*x)*x)));
  }
}

// ==================================================================================
void calcHessian(averages* path, parameters params){
  // fill out the Hessian variables for the structure
  // initial params must have the positions filled
  double x;
  int i;

  for(i=0;i<params.NumB;i++){
    x = path[i].pos;
    //path[i].Hessian= -(6.75 + x * (-10.125 + x * (-34.1719 + (56.9531 - 21.3574 * x) * x)));
    path[i].Hessian= -6.75 + x * (10.125 + x * (34.1719 + x * (-56.9531 + 21.3574 * x)));
  }
}

// ==================================================================================
void calcDeltae(averages* path, parameters params){
  // find the change in energy for each point along the path
  int i;
  for(i=0;i<params.NumB-1;i++){
    path[i].deltae = Pot(path[i+1].pos)-Pot(path[i].pos) + 0.5*(path[i+1].Force + path[i].Force)*(path[i+1].pos - path[i].pos) + params.deltat*0.25*( path[i+1].Force*path[i+1].Force - path[i].Force*path[i].Force);
  }
}


// ==================================================================================
void calcdg(averages* path, parameters params){
  // calculate dG/dx for the structure
  // intial positions/force/Hessian must be filled
  int i;

  double dPlus;
  double dMinus;
  for(i=1;i<params.NumB-1;i++){
    dPlus = (path[i+1].pos-path[i].pos)*params.invdt - path[i].Force;
    dMinus = (path[i].pos-path[i-1].pos)*params.invdt + path[i].Force;

    path[i].dg  = (2.0*path[i].Force - path[i-1].Force - path[i+1].Force)*params.invdt;
    path[i].dg += path[i].Hessian*(dPlus-dMinus);
    path[i].dg += copysign(1.0,path[i].deltae) * ((path[i].Force - path[i+1].Force)*params.invdt - path[i].Hessian*dPlus);
    path[i].dg += copysign(1.0,path[i-1].deltae) * ((path[i-1].Force - path[i].Force)*params.invdt - path[i].Hessian*dMinus);
    path[i].dg *=0.5;
  }

  // boundaries are set to zero
  //path[0].dg=0.0;
  //path[params.NumB-1].dg=0.0;

}

// ==================================================================================
void calcSPDErhs(averages* path, parameters params){
  // calculate the right hand side vector for the SPDE
  // including the boundary conditions from the LHS matrix multiplication
  // note that the i_th position is (i+1)
  // temp1: (I+rL)xOld
  // temp2: dg term
  // temp3: noise term
  int i;
  double temp1, temp2, temp3;

  for(i=1;i<params.NumB-1;i++){
    temp1 = params.r*path[i-1].pos + (1.-2.*params.r)*path[i].pos + params.r*path[i+1].pos;
    temp2 = params.deltatau*path[i].dg;
    temp3 = params.noisePref*path[i].randlist;
    path[i].rhs=temp1-temp2+temp3;
  }
  path[1].rhs+=params.r*path[0].pos;
  path[params.NumB-2].rhs+=params.r*path[params.NumB-1].pos;

}


// ==================================================================================
void calcMDrhs(averages* path0, averages* path1, parameters params){
  int i;
  
  for(i=1; i<params.NumB-1;i++){

    path1[i].rhs  = 2.0*(params.r*path1[i-1].pos + (1.0-2.0*params.r)*path1[i].pos + params.r*path1[i+1].pos);
    path1[i].rhs += -1.0*(-params.r*path0[i-1].pos + (1.0+2.0*params.r)*path0[i].pos - params.r*path0[i+1].pos);
    path1[i].rhs -= 2.0*params.deltatau*path1[i].dg;
  }

  path1[1].rhs += 2.0*params.r*path1[0].pos;
  path1[1].rhs += -params.r*path0[0].pos;
  path1[params.NumB-2].rhs += 2.0*params.r*path1[params.NumB-1].pos;
  path1[params.NumB-2].rhs += -params.r*path1[params.NumB-1].pos;

}

// ==================================================================================
void calcPhi(averages* path, parameters params){
  int i;

  double Psi;
  double F0;
  double F1;

  for(i=0;i<params.NumB-1;i++){
    F0=path[i].Force;
    F1=path[i+1].Force;
    Psi=0.25*(F1*F1)+0.25*(F0*F0)+0.5*(F1-F0)*params.invdt*(path[i+1].pos-path[i].pos);
    
    path[i].Phi=Psi+abs(path[i].deltae)*params.invdt;
  }
}

// ==================================================================================
double calcEnergyChange(averages* path0, averages* path1, parameters params){
  int i;
  //temp storage for the sum of Phi for path0
  double Phi0=0.0;
  //temp storage for the sum of Phi for path1
  double Phi1=0.0;
  // temp var for first term in notes
  double temp1=0.0;
  // temp var for second term in notes
  double temp2=0.0;
  // temp var for the inner term of the sum (Lx1-dg1+Lx0-dg0)
  double temp2inner=0.0;

  double invdt2=params.invdt*params.invdt;

  // calculate the sum of Phi for each path
  for(i=0;i<params.NumB-1;i++){
    Phi0 += path0[i].Phi;
    Phi1 += path1[i].Phi;
  }

  //calculate the first term in the notes and multiply by constant
  for(i=1;i<params.NumB-1;i++){
    temp1 += (path1[i].dg+path0[i].dg)*(path1[i].pos*path0[i].pos);
  }
  temp1 = 0.5*temp1;

  //calculate the second term in the notes, first calcuating the inner 
  //   term and then the entire term, and finally the constant multiply
  for(i=1;i<params.NumB-1;i++){
    //calculate the inner term of the sum
    temp2inner  = invdt2*(path1[i+1].pos-2.0*path1[i].pos+path1[i-1].pos);
    temp2inner += -path1[i].dg;
    temp2inner += invdt2*(path0[i+1].pos-2.0*path0[i].pos+path0[i-1].pos);
    temp2inner += -path0[i].dg;
    // incrimenting the sum
    temp2 += temp2inner*(path1[i].dg-path0[i].dg);
  }
  temp2 *= 0.25*params.deltatau;

  //return the change in energy for a single step
  return Phi1-Phi0-temp1-temp2;

}


// ==================================================================================
// Gaussian elimination for computing L.x=b where
//   L is tridiagonal matrix
//     mainDiag: 1+2r
//     upper(lower)Diag: -r
//
//   input: the filled path0 state 
//   output: new path1.pos  
void GaussElim(averages* path0, averages* path1, parameters params){

  int NumL = params.NumB-2;
  int i;
  double temp;
  
  // declare the L matrix arrays
  double al[NumL-1];
  double am[NumL];
  double au[NumL-1];

  //initialize the tridiagonal arrays
  for(i=0;i<NumL-1;i++){ 
    al[i]=-params.r;
    am[i]=(1.0+2.0*params.r);
    au[i]=-params.r;
  }
  am[NumL-1]=(1.0+2.0*params.r);


  // Gaussian Elimination
  for(i=0;i<NumL-1;i++){
    temp=-al[i]/am[i];
    am[i+1]+=au[i]*temp;
    path0[i+2].rhs+=path0[i+1].rhs*temp;
  }

  // Back substitution
  for(i=0;i<NumL-1;i++){
    temp=-au[NumL-i-2]/am[NumL-i-1];
    path0[NumL-i-1].rhs+=path0[NumL-i].rhs*temp;
  }
  
  // Divide by main diagonal
  for(i=0;i<NumL;i++){
    path0[i+1].rhs=path0[i+1].rhs/am[i];
  }

  // Now the rhs vector transformed to the new position vector
  // save the new positions to the path1 positions
  for(i=1;i<params.NumB-1;i++){
    path1[i].pos=path0[i].rhs;

  // set the boundary conditions on the new path positions
  // BCs are not included in the Gaussan elimination
  path1[0].pos=path0[0].pos;
  path1[params.NumB-1].pos=path0[params.NumB-1].pos;
  }
}

// ==================================================================================
void quadVar(averages* path, parameters params){
  // fill out the Hessian variables for the structure
  // initial params must have the positions filled
  double qv=0.0;
  int i;

  for(i=0;i<params.NumB-1;i++){
    qv+=(path[i+1].pos-path[i].pos)*(path[i+1].pos-path[i].pos);
  }
  printf("qv=%f\n",qv/(2.0*params.eps*params.deltat*(params.NumB-1)));
}
// ==================================================================================

