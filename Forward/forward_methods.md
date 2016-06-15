# Forward integration methods

This is the directory where all forward methods are stored.
By *forward method*, I am referring to any sampling code that makes a forward evolving trajectory rather than a double ended *path*.

All of these codes will be absorbed into the Ctypes code later and renamed into a single code.

## forward_Bs
Calculate the broad well fraction for poor mans doulble ended forward paths with boundary conditions.

## forward_int_stats
no boundary conditions but all integration methods.
used to calculate B(s) for {leapfrog,midpt,simpson} quadratures at different time steps.

## gen_path_broad_frac.py

supporting script used to create trajectories and only accept the path when B(s) is within a range and the ending position, x(T), is within a range.

## gen_path_ctypes

The code that will rule them all. 
Impliments ctypes to speed up the generation of trajectories with all three quadratures.
the ctypes part of the codde stops at the generation of the path and is taken over by normal python.
