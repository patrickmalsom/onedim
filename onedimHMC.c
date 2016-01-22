//// ==================================================================================
void calcPosBar(averages* path, parameters params){
  // fill out the Potential variables for the struce
  // initial params must have the positions filled
  int i;

  #pragma omp parallel for
  for(i=0;i<params.NumB-1;i++){
    path[i].posBar = (path[i].pos+path[i+1].pos)*0.5;
  }
}
// ==================================================================================
void calcPot(averages* path, parameters params){
  // fill out the Force variables for the struce
  // initial params must have the positions filled
  int i;
  #pragma omp parallel for
  for(i=0;i<params.NumB;i++){
    path[i].U = Pot(path[i].pos);
  }
}
// ==================================================================================
void calcForces(averages* path, parameters params){
  // fill out the Force variables for the struce
  // initial params must have the positions filled
  int i;
  #pragma omp parallel for
  for(i=0;i<params.NumB;i++){
    path[i].F = Force(path[i].pos);
  }
}
// ==================================================================================
void calcForcesBar(averages* path, parameters params){
  int i;
  #pragma omp parallel for
  for(i=0;i<params.NumB;i++){
    path[i].Fbar = Force(path[i].posBar);
  }
}
// ==================================================================================
void calcForcesPrime(averages* path, parameters params){
  int i;
  #pragma omp parallel for
  for(i=0;i<params.NumB;i++){
    path[i].Fp = ForcePrime(path[i].pos);
  }
}
// ==================================================================================
void calcForcesPrimeBar(averages* path, parameters params){
  int i;
  #pragma omp parallel for
  for(i=0;i<params.NumB;i++){
    path[i].Fpbar = ForcePrime(path[i].posBar);
  }
}

// ==================================================================================
void calcForcesDoublePrime(averages* path, parameters params){
  int i;
  #pragma omp parallel for
  for(i=0;i<params.NumB;i++){
    path[i].Fpp = ForceDoublePrime(path[i].pos);
  }
}
// ==================================================================================
void calcForcesDoublePrimeBar(averages* path, parameters params){
  int i;
  #pragma omp parallel for
  for(i=0;i<params.NumB;i++){
    path[i].Fppbar = ForceDoublePrime(path[i].posBar);
  }
}
// ==================================================================================
void calcG(averages* path, parameters params){
  int i;
  #pragma omp parallel for
  for(i=0;i<params.NumB;i++){
    path[i].G=0.5*path[i].F*path[i].F + params.eps*path[i].Fp;
  }
}
// ==================================================================================
void calcgradG(averages* path, parameters params){
  int i;
  #pragma omp parallel for
  for(i=0;i<params.NumB;i++){
    path[i].gradG=path[i].F*path[i].Fp + params.eps*path[i].Fpp;
  }
}

/* =====================================================================
 * Function: LInverse
 * -----------------------------
 * Calculates L^{-1} Grad G for a path struct
 *   - Struct potentials must be filled before calling (calcPotentials)
 *   - THIS LOOP IS RECURSIVE! take care when parallelizing
 * 
 *   path: input path_struct to calculate LinvG for
 *   params: (constant) parameters for the run
*/
void LInverse(averages* path, parameters params)
{
  double lasti;
  int n;
  int Numl=params.NumB-2;
  double vecdg[Numl];
  double veci0[Numl];
  double veci1[Numl];

  for(n=0;n<params.NumB-2;n++){
    vecdg[n]=path[n+1].gradG*params.deltat*params.deltat;
  }
  veci0[0]=vecdg[0];
  veci1[0]=vecdg[0];
  //this for loop is recursive!!!!
  for(n=1;n<Numl;n++){
    veci0[n]=veci0[n-1]+vecdg[n];
    veci1[n]=veci1[n-1]+((double)(n+1))*vecdg[n];
  }

  lasti=veci0[Numl-1]-(path[params.NumB-1].pos-path[0].pos+veci1[Numl-1])/((double)(Numl+1));

  for(n=0;n<Numl;n++){
    path[n+1].LinvG=path[0].pos+((double)(n+1))*(veci0[n]-lasti) - veci1[n];
  }
  path[0].LinvG=path[0].pos;
  path[params.NumB-1].LinvG=path[params.NumB-1].pos;

}

/* =====================================================================
 * Function: calcSPDEpos
 * -----------------------------
 * Calculate the new position for the MALA step (SPDE)
 *   - fills path1[i].pos with the new path positions
 *   - requires that path0 is completely initialized
 *     (calcPotentials, Liverse, BrownianBridge)
 * 
 *   path0: input path_struct (filled)
 *   path1: output path_struct
 *   params: (constant) parameters for the run
*/
void calcSPDEItopos(averages* path0, averages* path1, parameters params){
  // calculate the right hand side vector for the SPDE
  // this is x^{(1)} vector in the notes

  int i;
  double h=sqrt(2.0*params.deltatau);
  double si=(4.0*h)/(4.0+h*h);
  double co=(4.0-h*h)/(4.0+h*h);

  for(i=1;i<params.NumB-1;i++){
  path1[i].pos=si*h*0.5*path0[i].LinvG + si*path0[i].bb + co * path0[i].pos;
  }

  //set boundary conditions
  path1[0].pos=path0[0].pos;
  path1[params.NumB-1].pos=path0[params.NumB-1].pos;

}

// ==================================================================================
void calcSPDEFiniterhs(averages* path, parameters params){
  // calculate the right hand side vector for the SPDE
  // including the boundary conditions from the LHS matrix multiplication
  // note that the i_th position is (i+1)
  // temp1: (I+rL)xOld
  // temp2: dg term
  // temp3: noise term
  int i;
  double temp1, temp2, temp3;

  #pragma omp parallel for private(temp1,temp2,temp3)
  for(i=1;i<params.NumB-1;i++){
    temp1 = params.r*path[i-1].pos + (1.-2.*params.r)*path[i].pos + params.r*path[i+1].pos;
    temp2 = params.deltatau*path[i].dg;
    temp3 = params.noisePref*path[i].randlist;
    path[i].rhs=temp1-temp2+temp3;
  }
  path[1].rhs+=params.r*path[0].pos;
  path[params.NumB-2].rhs+=params.r*path[params.NumB-1].pos;

}

/* =====================================================================
 * Function: calcMDpos
 * -----------------------------
 * Calculate the new position for a molecular dynamics step (MD)
 *   - fills path2[i].pos with the new path positions
 *   - requires that path0 and path1 are completely initialized
 *     (calcPotentials, Liverse)
 * 
 *   path0: old path_struct (filled)
 *   path1: current path_struct (filled)
 *   path2: output path_struct
 *   params: (constant) parameters for the run
*/
void calcMDItopos(averages* path0, averages* path1, averages* path2, parameters params){
  // calculate the right hand side vector for the MD steps
  // this is x^{(2)} vector in the notes

  int i;
  double h=sqrt(2.0*params.deltatau);
  double si=(4.0*h)/(4.0+h*h);
  double co=(4.0-h*h)/(4.0+h*h);

  for(i=1;i<params.NumB-1;i++){
  path2[i].pos=(path1[i].pos-path0[i].pos) + (2.0*co - 1.0)*path1[i].pos + si*h*path1[i].LinvG;
  }

  //set boundary conditions
  path2[0].pos=path1[0].pos;
  path2[params.NumB-1].pos=path1[params.NumB-1].pos;

}

/* =====================================================================
 * Function: calcEnergyChange
 * -----------------------------
 * Calculate the energy error made when integrating from path0 to path1
 *   - requires that path0 and path1 are completely initialized
 *     (calcPotentials, Liverse)
 * 
 *   path0: current path_struct (filled)
 *   path1: new path_struct (filled)
 *   params: (constant) parameters for the run
 *
 *   return: the negative argument of the metropolis hastings test: dE/(2 epsilon)
*/
double calcEChangeIto(averages* path0, averages* path1, parameters params){

  int i;

  double h=sqrt(2.0*params.deltatau);
  double cot=(4.0l-h*h)/(4.0l*h);
  double csc=(4.0l+h*h)/(4.0l*h);

  double lambda1=0.0;
  double lambda2=0.0;
  // tempsum is a temporary double for storage of the integrals in the lambdas
  // it should be set to zero at the beginning of each integration evaluation.
  double tempsum;

  // first integral in Lambda1
  tempsum=0.0;
  for(i=0;i<params.NumB;i++){
    tempsum+= (path1[i].pos *path0[i].gradG - path0[i].pos*path1[i].gradG);
  }
  lambda1-= h/(4.0*params.eps)*csc * params.deltat*tempsum;

  // second integral in Lambda1
  tempsum=0.0;
  for(i=0;i<params.NumB;i++){
    tempsum+= (path0[i].pos*path0[i].gradG - path1[i].pos*path1[i].gradG);
  }
  lambda1+=h/(4.0*params.eps)*cot * params.deltat*tempsum;

  // third integral in Lambda1
  tempsum=0.0;
  for(i=0;i<params.NumB;i++){
    tempsum+= ( path0[i].gradG*path0[i].LinvG - path1[i].gradG*path1[i].LinvG );
  }
  lambda1+=h*h/(16.0*params.eps) * params.deltat*tempsum;

  // integral in Lambda2
  tempsum=0.0;
  for(i=0;i<params.NumB;i++){
    tempsum+= (path1[i].G - path0[i].G);
  }
  lambda2+=0.5/params.eps * params.deltat*tempsum;

  //return the change in energy for a path step
  return( lambda1+lambda2 );

}


/* =====================================================================
 * Function: generateBB
 * -----------------------------
 * Make a normalized Brownian bridge for use in the MALA (SPDE) step
 *   - saves the bridge to path.bb
 *   - requires that path.randlist is filled (from the python side)
 * 
 *   path0: any path_struct (randlist filled)
 *   params: (constant) parameters for the run
*/
void generateBB(averages* path, parameters params){
//generate a Brownian bridge and store to BrownianBridge
  int n;
  double xn;
  double doubleNUMu=((double)(params.NumB-1));

  double endPtCorr;
  double sum, term, term0;
  double alpha;

  double sqdu=sqrt(2.0*params.eps*params.deltat);
  path[0].bb=0.0l;
  for(n=1;n<params.NumB;n++)
  {
    path[n].bb=path[n-1].bb+sqdu*path[n-1].randlist;
  }

  xn=path[params.NumB-1].bb/((double)(params.NumB-1));
  for(n=1;n<params.NumB-1;n++)
  {
    path[n].bb-=((double)(n))*xn;
  }
  path[params.NumB-1].bb=0.0l;

  //normalize the brownian bridge to have accurate quad var
  endPtCorr= (path[0].bb-path[params.NumB-1].bb)*(path[0].bb-path[params.NumB-1].bb)/doubleNUMu;
  sum=0.0l;
  for(n=0;n<params.NumB-1;n++){
    sum+= (path[n].bb-path[n+1].bb)*(path[n].bb-path[n+1].bb);
  }
  alpha=sqrt((doubleNUMu*params.deltat*2.0*params.eps-endPtCorr)/(sum-endPtCorr));
  term=(1.0l-alpha)*(path[params.NumB-1].bb-path[0].bb)/doubleNUMu;
  term0=(1.0l-alpha)*path[0].bb;
  //term and term0 have a subtraction of roughly equal numbers and thus is not very accurate
  // alpha is ~1 with an error of 10^-4 or 5 for sample configs. This makes the routine
  //nondeterministic between Fortran and C
  for(n=1;n<params.NumB-1;n++){
    path[n].bb=alpha*path[n].bb+term0+((double)(n-1))*term;
  }
}

// ==================================================================================
void calcDeltae(averages* path, parameters params){
  // find the change in energy for each point along the path
  int i;

  //midpt (finite method -> 1)
  if(params.method == 1){
    #pragma omp parallel for
    for(i=0;i<params.NumB-1;i++){
      path[i].deltae = path[i+1].U - path[i].U + path[i].Fbar*(path[i+1].pos -path[i].pos);
    }
  }

  //leapfrog (finite method -> 2)
  else if(params.method == 2){
    #pragma omp parallel for
    for(i=0;i<params.NumB-1;i++){
      path[i].deltae = path[i+1].U - path[i].U + (path[i+1].pos-path[i].pos) * (path[i+1].F+path[i].F)*0.5 + params.deltat*0.25*(path[i+1].F*path[i+1].F - path[i].F*path[i].F);
    }
  }

  //simpsons (finite method -> 3)
  if(params.method == 3){
    #pragma omp parallel for
    for(i=0;i<params.NumB-1;i++){
      path[i].deltae = path[i+1].U - path[i].U + (path[i].F+4.0*path[i].Fbar+path[i+1].F)*(path[i+1].pos -path[i].pos)*0.1666666666666666666;
    }
  }
}

// ==================================================================================
void calcMDFiniterhs(averages* path0, averages* path1, parameters params){
  int i;
  
  #pragma omp parallel for
  for(i=1; i<params.NumB-1;i++){

    path1[i].rhs  = 2.0*(params.r*path1[i-1].pos + (1.0-2.0*params.r)*path1[i].pos + params.r*path1[i+1].pos);
    path1[i].rhs += -1.0*(-params.r*path0[i-1].pos + (1.0+2.0*params.r)*path0[i].pos - params.r*path0[i+1].pos);
    path1[i].rhs -= 2.0*params.deltatau*path1[i].dg;
  }

  // remember that the path1[0].pos is the boundary condition of the path2 struct
  // even though it is unknown at this point, the BC is known
  path1[1].rhs += params.r*path0[0].pos;
  path1[params.NumB-2].rhs += params.r*path1[params.NumB-1].pos;

}


// ==================================================================================
double calcEnergyChangeFinite(averages* path0, averages* path1, parameters params){
// thesis_eqn 7.54 
  int i;
  //temp storage for the sum of Phi for path0
  double Phi0=0.0;
  //temp storage for the sum of Phi for path1
  double Phi1=0.0;

  // temp var for first sum (dg1+dg0)*(x1-x0)
  double temp1=0.0;

  // temp var for second sum
  double temp2=0.0;
  // temp var for the first inner term of the second sum (Lx1-dg1+Lx0-dg0)
  double temp2inner=0.0;

  double invdt2=params.invdt*params.invdt;

  // calculate the sum of Phi for each path
  for(i=0;i<params.NumB-1;i++){
    Phi0 += path0[i].Phi;
    Phi1 += path1[i].Phi;
  }

  //calculate the first sum in thesis_eqn 7.54 and multiply by half
  for(i=1;i<params.NumB-1;i++){
    temp1 += (path1[i].dg+path0[i].dg)*(path1[i].pos-path0[i].pos);
  }
  temp1 = 0.5*temp1;

  //calculate the second term in the notes, first calcuating the inner 
  //   term and then the entire term, and finally the constant multiply
  for(i=1;i<params.NumB-1;i++){
    //calculate the inner term of the sum
    temp2inner  = invdt2*(path1[i+1].pos-2.0*path1[i].pos+path1[i-1].pos);
    temp2inner += invdt2*(path0[i+1].pos-2.0*path0[i].pos+path0[i-1].pos);
    temp2inner += -path1[i].dg;
    temp2inner += -path0[i].dg;
    // incrimenting the sum
    temp2 += temp2inner*(path1[i].dg-path0[i].dg);
  }
  temp2 *= 0.25*params.deltatau;

  //return the change in energy for a single step
  // this returns the negative of the entire MHMC test argument (dE*dt/2/eps)
  return (Phi1-Phi0-temp1-temp2)*params.deltat/(2.0*params.eps);

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
// calcdg: Calculate the sum of the energy error between at each time interval
// (for a specific itegration method) for the entire path
void calcdg(averages* path, parameters params){
  // calculate dG/dx for the structure
  // intial positions/force/Hessian must be filled
  int i;
  double sixth=0.1666666666666666666;
  double third=0.3333333333333333333;
  double Jn; //J_n for simpsons rule
  double Jnm1; // J_{n-1} for simpsons rule
  double Fwn; // weighted force F_w(x_n)
  double Fwnm1; // weighted force F_w(x_{n-1})
  double dFwn; // deriv of weighted force F_w(x_n)
  double dFwnm1; // deriv of weighted force F_w(x_{n-1})

  //midpt (method -> 1)
  if(params.method == 1){
    #pragma omp parallel for
    for(i=1;i<params.NumB-1;i++){
      path[i].dg  = 0.5 *(path[i].Fbar *path[i].Fpbar + params.eps * path[i].Fppbar /(1.0 - 0.5* params.deltat * path[i].Fpbar));
      path[i].dg += 0.5 *(path[i-1].Fbar *path[i-1].Fpbar + params.eps * path[i-1].Fppbar /(1.0 - 0.5* params.deltat * path[i-1].Fpbar));
      path[i].dg -= params.invdt * ( path[i].F - path[i].Fbar + 0.5*(path[i+1].pos - path[i].pos) * path[i].Fpbar );
      path[i].dg -= params.invdt * ( -path[i].F + path[i-1].Fbar + 0.5*(path[i].pos - path[i-1].pos) * path[i-1].Fpbar );
    }
  }

  //leapfrog (method -> 2)
  else if(params.method == 2){
    #pragma omp parallel for
    for(i=1;i<params.NumB-1;i++){
      //direct calculation
      path[i].dg = path[i].F*path[i].Fp + params.invdt*(path[i].F - path[i-1].F - path[i].Fp*(path[i+1].pos - path[i].pos));
      //thesis calculation
      //path[i].dg  = path[i].Fp*path[i].F - 0.5*params.invdt*(path[i+1].F - path[i].F) ;
      //path[i].dg += 0.5*params.invdt*(path[i].F - path[i-1].F);
      //path[i].dg -= params.invdt * ( path[i].F - 0.5*(path[i+1].F + path[i].F) + 0.5*(path[i+1].pos - path[i].pos)*path[i].Fp);
      //path[i].dg -= params.invdt * (-path[i].F + 0.5*(path[i].F + path[i-1].F) + 0.5*(path[i].pos - path[i-1].pos)*path[i].Fp);
    }
  }

  //simpsons (finite method -> 3)
  if(params.method == 3){
    #pragma omp parallel for
    for(i=1;i<params.NumB-1;i++){

      // calculate the Jacobian terms
      Jn = 1.0 - params.deltat*sixth*(path[i+1].Fp + 2.0*path[i].Fpbar);
      Jnm1 = 1.0 - params.deltat*sixth*(path[i].Fp + 2.0*path[i-1].Fpbar);
      Fwn = sixth*(path[i].F + 4.0*path[i].Fbar + path[i+1].F);
      Fwnm1 = sixth*(path[i-1].F + 4.0*path[i-1].Fbar + path[i].F);
      dFwn = sixth*(path[i].Fp+2.0*path[i].Fpbar);
      dFwnm1 = sixth*(path[i].Fp+2.0*path[i-1].Fpbar);

      path[i].dg = dFwn*(Fwn- params.invdt*(path[i+1].pos - path[i].pos)) + third*params.eps*path[i].Fppbar/Jn;
      path[i].dg+= dFwnm1*(Fwnm1- params.invdt*(path[i].pos - path[i-1].pos)) +  + third*params.eps*(path[i].Fpp + path[i-1].Fppbar)/Jnm1;
      path[i].dg +=params.invdt*( Fwn- Fwnm1);
      //path[i].dg-= params.invdt * (path[i].F + Fwn + sixth * (path[i+1].pos - path[i].pos) * (path[i].Fp + 2.0 * path[i].Fpbar));
      //path[i].dg-= params.invdt * (-path[i].F + Fwnm1 + sixth * (path[i].pos - path[i-1].pos) * (path[i].Fp + 2.0 * path[i-1].Fpbar));
    }
  }
}

// ==================================================================================
//defined from eqn 7.27 and 7.28 as well
void calcPhi(averages* path, parameters params){
  int i;
  double Fw;
  double sixth=0.1666666666666666666;

  //midpt (method -> 1)
  if(params.method == 1){
    #pragma omp parallel for
    for(i=0;i<params.NumB-1;i++){
      path[i].Phi = 0.5*path[i].Fbar*path[i].Fbar - 2.0 * params.eps * params.invdt * log(1.0 - params.deltat*0.5*path[i].Fpbar) - path[i].deltae*params.invdt;
    }
  }

  //leapfrog (method -> 2)
  else if(params.method == 2){
    #pragma omp parallel for
    for(i=0;i<params.NumB-1;i++){
      //using deltae[i]
      //path[i].Phi = 0.25*path[i+1].F*path[i+1].F + 0.25*path[i].F*path[i].F + 0.5*(path[i+1].pos - path[i].pos)*params.invdt * (path[i+1].F - path[i].F) - path[i].deltae*params.invdt;
      //direct calculation
      path[i].Phi = 0.5*path[i].F*path[i].F - path[i].F*(path[i+1].pos - path[i].pos)*params.invdt;
    }
  }

  //simpsons (finite method -> 3)
  if(params.method == 3){
    #pragma omp parallel for
    for(i=0;i<params.NumB-1;i++){
      Fw=sixth*(path[i].F+4.0*path[i].Fbar+path[i+1].F);
      //path[i].Phi = 0.5*Fw*Fw - 2.0 * params.eps * params.invdt * log(1.0 - params.deltat*sixth*(path[i+1].Fp + 2.0*path[i].Fpbar));
      path[i].Phi = Fw*(0.5* Fw - (path[i+1].pos - path[i].pos)*params.invdt) - 2.0 * params.eps * params.invdt * log(1.0 - params.deltat*sixth*(path[i+1].Fp + 2.0*path[i].Fpbar));

    }
  }
}

/* =====================================================================
 * Function: quadVar
 * -----------------------------
 * Calculate the quadratic variation of a path
 * 
 *   path: any path_struct (pos filled)
 *   params: (constant) parameters for the run
 *
 *   return: normalized quadratic variation 
 *           sum(x_1-x_0)^2 /(2 eps T) -> 1.0
*/
double quadVar(averages* path, parameters params){

  double qv=0.0;
  int i;

  for(i=0;i<params.NumB-1;i++){
    qv+=(path[i+1].pos-path[i].pos)*(path[i+1].pos-path[i].pos);
  }
  return(qv/(2.0*params.eps*params.deltat*(params.NumB-1)));
}

