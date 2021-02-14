# The MPLAPACK is a multiprecision linear algebra package based on BLAS and LAPACK.
This package is rewritten in C++, and supports several high precision
libraries like GMP, MPFR and QD etc so that users can choose for user's
convenience. The MPLAPACK is a free software (2-clause BSD style license with
original license by LAPACK).

# How to build and install
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack ; autoreconf --force --install ; aclocal ; autoconf ; automake; autoreconf --force --install
$ ./configure --enable-gmp --enable-__float128 --enable-dd --enable-mpfr
```
do not worry about a lot of warnings in the second commands.

# Docker build

```
$ docker image build -t mplapack:latest -f Dockerfile_ubuntu20.04 .
$ docker run -it mplapack:latest
```

# Docker build + fable (fortran to c converter)
https://github.com/cctbx/cctbx_project/tree/master/fable
FABLE: Automatic Fortran to C++ conversion (https://doi.org/10.1186/1751-0473-7-5).

## How to build
```
docker image build -t mplapack:f2c_fable -f Dockerfile_fable .
```

## How to convert a sample fortran to C++
```
$ docker run -it mplapack:f2c_fable
$ cd ; fable.cout sample.f
```

# Docker build for rename2 mplapack branch

## How to build
```
docker image build -t mplapack:rename2mplapack -f Dockerfile_rename2mplapack .
```

# Acknowledgement:

This work has been supported by:
The Special Postdoctoral Researchers' Program of RIKEN (2008, 2009)
Grant-in-Aid for Scientific Research (B) 21300017 from the Japan Society for the Promotion of Science (2009, 2010, 2011).
Microsoft Research CORE6 (2010). 

# Special thanks to:

Fujisawa, Katsuki
Goto, Kazushige
Himeno, Ryutaro
.

Please refer userman.pdf for details. Please enjoy!

# contact
NAKATA Maho <maho@riken.jp>
