/*   ItoHMC.c 
     Written Winter 2014 -- Patrick Malsom
 
C library for the Ito HMC algorithm
Functions used to generate the SPDE and the MD steps

This is a shared library that is linked to python (ItoHMC.py)
*/

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
typedef struct _path
{
  // positions
  double pos;
  double posBar;
  // random gaussian numbers
  double randlist;
  // forces
  double F;
  double Fbar;
  double Fp;
  double Fpbar;
  double Fpp;
  double Fppbar;
  // Finite
  double deltae;
  double dg;
  double Phi;
  double rhs;
  // Ito
  double bb;
  double G;
  double gradG;
  double LinvG;
} averages;

