# The MPLAPACK is a multiprecision linear algebra package based on BLAS and LAPACK.
This package is rewritten in C++, and supports several high precision
libraries like GMP, MPFR and QD etc so that users can choose for user's
convenience. The MPLAPACK is a free software (2-clause BSD style license with
original license by LAPACK).

# Capabilities
* MPBLAS: All BLAS routines can be done in multiple precision arithmetic.
* MPLAPACK: in ver 1.0.0 all Real version of following solvers (complex version will be added in 2.0)
* Linear Equations
* Linear Least Squares (LLS) Problems
* Generalized Linear Least Squares (LSE and GLM) Problems
* Standard Eigenvalue and Singular Value Problems
* Symmetric Eigenproblems (SEP)
* Nonsymmetric Eigenproblems (NEP)
* Singular Value Decomposition (SVD) 
* Generalized Eigenvalue and Singular Value Problems
* Generalized Symmetric Definite Eigenproblems (GSEP)
* Generalized Nonsymmetric Eigenproblems (GNEP)
* Generalized Singular Value Decomposition (GSVD) 

# Supported multiple precision libraries and floating point formats

* MPFR + MPC https://www.mpfr.org/ and http://www.multiprecision.org/mpc/  (arbitrary precision with IEEE like rounding mode)
* GMP https://gmplib.org/ (arbitrary precision)
* double (binary64)
* DD, QD (https://www.davidhbailey.com/dhbsoftware/) (DD=approx. binary128, QD=approx. binary256)
* _Float128 (binary128; via glibc or libquadmath; automatically detected)
* _Float64x (extended precision of double; binary80 Intel CPU only)

We use MPFR + MPC as the primary arithmetic class.

# Supported platforms

* Ubuntu 20.04, 18.04 (amd64, AArch64)
* CentOS 7,8 (amd64, AArch64)
* Ubuntu 20.04 (amd64) + Intel oneAPI
* macOS (Intel) + macports (you may use homebrew instead, small modification of build script req'ed)
* Windows (64bit; mingw64 on Ubuntu with wine64)

# How to build on Linux and Win (using Docker; recommended)

Ubuntu 20.04 (amd64, AArch64)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ /usr/bin/time docker build -t mplapack:ubuntu2004 -f Dockerfile_ubuntu20.04 . 2>&1 | tee log.ubuntu2004
```

CentOS 7 (amd64, AArch64)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ /usr/bin/time docker build -t mplapack:centos7 -f Dockerfile_CentOS7 . 2>&1 | tee log.CentOS7
```

CentOS 8 (amd64, AArch64)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ /usr/bin/time docker build -t mplapack:centos8 -f Dockerfile_CentOS8 . 2>&1 | tee log.CentOS8
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
$ git clone https://github.com/nakatamaho/mplapack.git -b v0.9.3 --depth 1
$ cd mplapack
$ pushd mplapack/debug ; bash gen.Makefile.am.sh ; popd
$ autoreconf --force --install ; aclocal ; autoconf ; automake; autoreconf --force --install
$ CXX="g++-mp-9" ; export CXX
$ CC="gcc-mp-9" ; export CC
$ FC="gfortran-mp-9"; export FC
$ F77="gfortran-mp-9"; export F77
$ ./configure --prefix=/Users/maho/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-debug=yes
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

## MPLAPACK 1.0.0 Release Process 
This is the release schedule for MPLAPACK 1.0.0

| Action | Date | Status | Description |
| --- | --- | --- | --- |
| QA of blas (mpfr)                    | 2021-08-09 |100% | compare with original BLAS | 
| QA of blas (gmp)                     | 2021-08-09 |100% | compare with MPFR MPBLAS   | 
| QA of blas (_Float128)               | 2021-08-09 |100% | compare with MPFR MPBLAS   | 
| QA of blas (_Float64x)               | 2021-08-09 |100% | compare with MPFR MPBLAS   | 
| QA of blas (double)                  | 2021-08-09 |100% | compare with MPFR MPBLAS   | 
| QA of blas (dd)                      | 2021-08-09 |100% | compare with MPFR MPBLAS   | 
| QA of blas (qd)                      | 2021-08-09 |100% | compare with MPFR MPBLAS   | 
| QA of lapack (mpfr)                  | 2021-08-15 |100% | compare with original LAPACK | 
| QA of lapack (gmp)                   | 2021-08-20 |100% | compare with MPFR LAPACK |  
| QA of lapack (_Float128)             | 2021-08-15 |100% | compare with MPFR LAPACK |  
| QA of lapack (_Float64x)             | 2021-08-19 |100% | compare with MPFR LAPACK |  
| QA of lapack (double)                | 2021-08-16 |100% | compare with MPFR LAPACK |  
| QA of lapack (dd)                    | 2021-08-19 |100% | compare with MPFR LAPACK |  
| QA of lapack (qd)                    | 2021-08-16 |100% | compare with MPFR LAPACK |  
| QA of lin (mpfr, real)               | 2021-08-15 |100% |                          | 
| QA of lin (gmp, real)                | 2021-08-20 |100% |                          | 
| QA of lin (_Float128, real)          | 2021-08-15 |100% |                          | 
| QA of lin (_Float64x, real)          | 2021-08-15 |100% |                          | 
| QA of lin (double, real)             | 2021-08-15 |100% |                          | 
| QA of lin (dd, real)                 | 2021-08-15 |100% |                          | 
| QA of lin (qd, real)                 | 2021-08-17 |100% |                          | 
| QA of eig (mpfr, real)               | 2021-08-18 |100% |                          | 
| QA of eig (gmp, real)                | 2021-08-20 |100% | errors of eig/Red.in are bit large but ignorable |
| QA of eig (_Float128, real)          | 2021-08-17 |100% |                          | 
| QA of eig (_Float64x, real)          | 2021-08-17 |100% |                          | 
| QA of eig (double, real)             | 2021-08-17 |100% |                          | 
| QA of eig (dd, real)                 | 2021-08-21 |100% |                          |
| QA of eig (qd, real)                 | 2021-08-21 |90%  |svd.in Underflow occurs   | 
| Build on Ubuntu 20.04 amd64          | 2021-08-13 |3e51d57| _Float128: supp. by libc only binary128, _Float64x == long double |
| Build on Ubuntu 18.04 amd64          | 2021-08-13 |3e51d57| _Float128: supp. by libc only binary128, _Float64x == long double |
| Build on Ubuntu 20.04 Intel oneAPI   | 2021-08-13 |3e51d57| _Float128: supp. by libc only binary128, _Float64x == long double | 
| Build on Ubuntu 20.04 mingw64        | 2021-08-13 |3e51d57| _Float128: supp. by libquadmath, _Float64x == __float80 == long double |
| Build on CentOS7 amd64               | 2021-08-13 |3e51d57| _Float128: supp. by libquadmath, _Float64x == __float80 == long double |
| Build on CentOS8 amd64               | 2021-08-13 |3e51d57| _Float128: supp. by libc only binary128, _Float64x == long double | 
| Build on MacOS Big Sur amd64         | 2021-08-13 |3e51d57| _Float128: supp. by libquadmath, _Float64x == __float80 == long double | 
| Build on Ubuntu 20.04 AArch64        | 2021-08-13 |3e51d57| _Float128: supp. by libc and _Float128 == long double |
| Build on CentOS7 AArch64             | 2021-08-13 |3e51d57| long double == binary128, _Float128 nor __float128 are supported |
| Build on CentOS8 AArch64             | 2021-08-13 |3e51d57| _Float128: supp. by libc and _Float128 == long double |
| Build on Ubuntu 20.04 PPC64LE        | 2021-08-13 |3e51d57| _Float128: supp. by libc and _Float128 != long double | 
| Add documents                        |            |     |                          | 

## MPLAPACK 2.0.0 Release Process 
This is the release schedule for MPLAPACK 2.0.0
| Action | Date | Status | Description |
| --- | --- | --- | --- |
| Impliment RFP version                |            |     |                          | 
| Impliment complex lin                |            |     |                          | 
| Impliment complex eig                |            |     |                          | 
| Impliment mixed precision version    |            |     |                          | 
| Fix Rlamch for GMP                   |            |     |                          | 
| Fix Rlamch for QD                    |            |     | svd.in fails.used with caution | 
| cleanup pow (REAL, long int)         |            |     |                          | 
| Get rid of compiler warnings         |            |     |                          | 
| Add more examples                    |            |     |                          | 

## MPLAPACK 3.0.0 Release Process 
This is the release schedule for MPLAPACK 3.0.0
| Action | Date | Status | Description |
| --- | --- | --- | --- |
| Impliment faster MPFR C++ wrapper like gmpxx.h |            |     |                          | 
| Drop GMP version                               |            |     | Since trigonometric functions req'ed | 
| optimized implimentations                      |            |     |                          | 

# History
* 2021/4/11 0.9.3 release. CentOS7 AArch64 support
* 2021/4/6  0.9.1 release. CentOS support
* 2021/4/1  0.9.0 release. Rename to mplapack. You must rename include files, etc. Rewrite and recompilation required.
* 2012/12/25: MPACK 0.8.0 NVIDIA C2050 support for Rgemm in double-double, and preliminary support for Intel Xeon Phi. ANNOUNCE
* 2012/12/20: MPACK 0.8.0-RC2 Build fixes on various platforms.
* 2012/12/05: Our Rgemm dd paper "A Fast implementation of matrix-matrix product in double-double precision on NVIDIA C2050 and application to semidefinite programming" is selected as the Best Papers of The Third International Conference on Networking and Computing. Slide is here.
* 2012/11/29: MPACK 0.8.0-RC1 CUDA version of Rgemm in double-double precision is integrated.
* 2012/10/13: CUDA 4.2 or later version of accelarated Rgemm in double-double precision on NVIDIA C2050 GPU is now available. Note it does not work on CUDA 4.1. Origial release announce is here., and preprint is available from here, and it will be presented at The Third International Conference on Networking and Computing Okinawa, Japan, December 5-7, 2012 .
* 2012/06/16: MPACK 0.7.0! Announcement
* 2012/06/16: Development has been migrated to SVN repository.
* 2011/10/28: Rgemm accleration in double-double precision on NVIDIA C2050 GPU is now available. Even though papers are not published, you can just try by "make". Note that only CUDA 3.2 is supported. Origial release announce is here.
* 2011/08/24: Rgemm accleration on NVIDIA C2050 GPU is coming. Unforutnately paper are rejected, so please wait... Here is a pdf slide.
* 2010/08/20: MPACK 0.6.7! Includes condition number estimators; Rgecon and Rpocon. Now 100 MLAPACK routines, and license has been changed to 2-caluse BSD style license. No longer use LGPLv3.
* 2010/08/6: MPACK 0.6.6! Build fix release. Tested on various Linux distributions.
* 2010/05/31: A paper for MPACK (0.6.4) in Japanese has been uploaded.
* 2010/05/21: MPACK 0.6.5! MPFR support, and MBLAS license has been changed to BSD style. Still MLAPACK part is LGPLv3. I'll replace hopefully soon.
* 2010/01/13: MPACK 0.6.4! BUILD FIX RELEASE! PLEASE CHECK ON YOUR ENVIRONMET! THANKS! ALSO WINDOWS IS SUPPORTED!
* 2009/12/18: POSTER and SLIDE ARE UPLOADED; MPACK (MBLAS/MLAPACK) poster at HPCS2010 in English, and slide in Japanese, and I did two seminars about MPACK (MBLAS/MLAPACK) @NII and @Tokyo Univ; here is the slide.
* 2009/11/24: MPACK 0.6.0!
* 2009/11/7: Add example page.
* 2009/10/9: MPACK 0.5.2, build fix on Intel Mac.
* 2009/10/6: The cvs version has faster Raxpy using OpenMP parallelism.
* 2009/10/5: MPACK 0.5.1, just minor updates.
* 2009/9/24: MPACK 0.5.0!
* 2009/9/17: GMP/QD/DD integration is going very well; now builds the mlapack part as well.
* 2009/9/11: Now CVS version of MBLAS supports QD and DD. I abandoned "explicit generation of instances" for bad performance.
* 2009/5/25: Now I switched to template programming so that the library will be most portable. We automatically generate GMP, QD and DD versions via "explicit generation of instances". For real calculations, source code level compatibility is retained and also we can made some optimized version as well, whereas complex code we need a small changes. Typically 2-3 percent performance loss have been observed.
* 2009/3/21: QD and DD version of SDPA (SDPA-QD, DD) have been uploaded. These packages include some part of MPACK in QD and DD.
* 2009/2/10: mpack-devel ML has been launched.
* 2009/2/10: mpack-0.0.9.tar.gz. 
* 2009/2/5: SDPA-GMP 7.1.2 has been released now is supported by MPACK (MBLAS/MLAPACK!)
* 2009/1/8: mpack-0.0.8.tar.gz. Moved to souceforge.net.
* 2008/8/22: mpack-0.0.6.tar.gz. Build fix release.
* 2008/8/14: mpack-0.0.5.tar.gz. Eighteen LAPACK routines (Rlamch, Rlae2, Rlaev2, Rlassq, Rlanst, Rlansy, Rlapy2, Rlarf, Rlarfg, Rpotf2, Rlarfb, Rlaset, Rlarft, Rlartg, Rlascl, Rlasrt, Rlasr, Rlatrd) with check programs and complete "blas.h" and "lapack.h" (prototypes for the BLAS and LAPACK) have been made.
* 2008/7/23: mpack-0.0.4.tar.gz. All test cases are imported.
* 2008/7/18: mpack-0.0.3.tar.gz. Fix copyright issues, make check partially works, except for MacOSX.
* 2008/7/17: mpack-0.0.2.tar.gz. Installation instructions on Fedora9 and Ubuntu 8.04 are now available. Build fix with gcc-4.3.
* 2008/7/15: mpack-0.0.1.tar.gz. Now configurable and installable :)
* 2008/7/11: Change GPL to LGPLv3 and upload mpack-0.0.0.tar.gz.
* 2008/6/27: Add mBSD header of original authors.
* 2008/6/26: Tar ball has been provided. All sources for MBLAS have been uploaded.
* 2008/6/24: This page has been created. 

# Oldpage
http://mplapack.sourceforge.net/

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

