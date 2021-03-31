# The MPLAPACK is a multiprecision linear algebra package based on BLAS and LAPACK.
This package is rewritten in C++, and supports several high precision
libraries like GMP, MPFR and QD etc so that users can choose for user's
convenience. The MPLAPACK is a free software (2-clause BSD style license with
original license by LAPACK).

# Supported platforms

* Ubuntu 20.04, 18.04 (amd64, AArch64)
* Ubuntu 20.04 (amd64) + Intel oneAPI
* macOS (Intel) + macports
* Windows (64bit)

# Supported multiple precision libraries and floating point formats

* MPFR + MPC https://www.mpfr.org/ and http://www.multiprecision.org/mpc/
* GMP https://gmplib.org/
* double (binary64)
* DD, QD (https://www.davidhbailey.com/dhbsoftware/)
* _Float128 (binary128; via glibc or libquadmath; automatically detected)
* _Float64x (extended precision of double; Intel CPU only)

We use MPFR + MPC as the primary arithmetic class.

# How to build on Linux and Win (using Docker; recommended)

Ubuntu 20.04 (amd64, AArch64)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ /usr/bin/time docker build -t mplapack:ubuntu2004 -f Dockerfile_ubuntu20.04 . 2>&1 | tee log.ubuntu2004
```

Ubuntu 18.04 (amd64, AArch64)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ /usr/bin/time docker build -t mplapack:ubuntu1804 -f Dockerfile_ubuntu18.04 . 2>&1 | tee log.ubuntu1804
```

Ubuntu 20.04 (using Intel oneAPI)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ /usr/bin/time docker build -t mplapack:ubuntu2004intel -f Dockerfile_ubuntu20.04_intel . 2>&1 | tee log.ubuntu2004.intel

```

Windows 64bit (using cross compiler on Ubuntu) 
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ /usr/bin/time docker build -t mplapack:mingw -f  Dockerfile_ubuntu20.04_mingw64 . 2>&1 | tee log.mingw
```

# How to build and install (on Mac, using macports)

```
$ sudo port install gcc9 coreutils git ccache
$ git clone https://github.com/nakatamaho/mplapack.git -b v0.9.0 --depth 1
$ cd mplapack
$ pushd mplapack/debug ; bash gen.Makefile.am.sh ; popd
$ autoreconf --force --install ; aclocal ; autoconf ; automake; autoreconf --force --install
$ ./configure --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-debug=yes
```

Note: Float64x is supported only on Intel CPUs.

# Docker build + FABLE (Automatic Fortran to C++ converter)

https://github.com/cctbx/cctbx_project/tree/master/fable
FABLE: Automatic Fortran to C++ conversion (https://doi.org/10.1186/1751-0473-7-5).

## How to build
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ /usr/bin/time docker image build -t mplapack:fable -f Dockerfile_fable . 2>&1 | tee log.fable
```

## How to convert a sample Fortran to C++
```
$ docker run -it mplapack:fable
$ cd ; fable.cout sample.f
```

# History
## 2021/4/1  0.9.0 release. Rename to mplapack. You must rename include files, etc. Recompilation required.

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
