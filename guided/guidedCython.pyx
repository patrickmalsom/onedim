from libc.math cimport exp, tanh, sqrt

def HH(double t, int n, double tH, double gamma):
    # normPref is 1/Sqrt[2^n n! Sqrt[Pi]]
    cdef double* normPref = [0.751125544464942, 0.531125966013598, 0.265562983006799, 0.108415633823010, 0.0383307149314439, 0.0121212363525988]
    cdef double* HHPoly = [ 1. , 0.  , 0.  , 0.   , 0. , 0. , \
                            0. , 2.  , 0.  , 0.   , 0. , 0. , \
                            -2., 0.  , 4.  , 0.   , 0. , 0. , \
                            0. , -12., 0.  , 8.   , 0. , 0. , \
                            12., 0.  , -48., 0.   , 16., 0. , \
                            0. , 120., 0.  , -160., 0. , 32. ]
    cdef double tempSum = 0.0
    for k in range(6):
        tempSum +=  sqrt(gamma) * (t-tH)**k * HHPoly[(6*n) + k]
    return gamma**(0.25) * normPref[n] * tempSum
