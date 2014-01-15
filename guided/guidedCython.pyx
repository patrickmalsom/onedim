from libc.math cimport exp, tanh

def m(double t):
    return 0.969624*tanh(2.0 * (t-4.8))

def A(double t):
    return -7.0*exp( -4.0 * (t - 5.0) * (t-5.0) )+7.52137

def U0(double x, double mean, double width):
    return 0.5*( (x - mean) * (x - mean) ) * width

def U(double x):
    return (x*x - 1.0)**2.0

def DeltaU(double x, double mean, double width):
    return U(x) - U0(x,mean,width)

def genStep(double x0, double v0h, double deltat, double thalf):
    return (x0 + v0h - 0.5*deltat*x0*A(thalf) + m(thalf)*A(thalf)*deltat) / (1. + 0.5*deltat*A(thalf))

def energyChange(double x0, double x1, double t0, double t1, double thalf):
    return U0(x1,m(t1),A(t1)) - U0(x0,m(t0),A(t0)) - (x1-x0)*A(thalf)*(x1+x0-2.0*m(thalf))*0.5
