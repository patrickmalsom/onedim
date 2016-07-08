# GSL

This simulation code uses gsl to generate random numbers.
Because gsl is not available on all systems, the code is statically compiled.
A statically compiled glibc math library (libm.a) is also required, and provided on fedora/centos via `glibc-static` package.

This code is *NOT* included in the repository and the build will fail without it!

GSL library (used v1.16 for published results)
[http://ftp.wayne.edu/gnu/gsl/]

Move to the forward_B project directory (where the makefile is located)

```
# get the gsl-1.16 code 
curl -O http://ftp.wayne.edu/gnu/gsl/gsl-1.16.tar.gz
# extract the tarball
tar -xvzf gsl-1.16.tar.gz
# move to the gsl dir
cd gsl-1.16
# configure the make with no shared libraries and build dir is current dir
./configure --disable-shared --prefix $(pwd)
# compile the code
make
# install in the current dir (can be changed with above --prefix)
make install
```

now go down a directory to the forward_B code and `make`
