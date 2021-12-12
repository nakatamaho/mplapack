# The MPLAPACK is a multiprecision linear algebra package based on BLAS and LAPACK.
This package is rewritten in C++, and supports several high precision
libraries like GMP, MPFR and QD etc so that users can choose for user's
convenience. The MPLAPACK is a free software (2-clause BSD style license with
original license by LAPACK).

# Capabilities
* MPBLAS: All BLAS routines can be done in multiple precision arithmetic.
* MPLAPACK: in ver 1.0.1 all Real version of following solvers (complex version will be added in 2.0)
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

# Manual
* https://arxiv.org/abs/2109.13406
* https://raw.githubusercontent.com/nakatamaho/mplapack/master/doc/manual/manual.pdf (updated frequently)
```
@misc{2109.13406,
Author = {Maho Nakata},
Title = {MPLAPACK version 1.0.0 user manual},
Year = {2021},
Eprint = {arXiv:2109.13406},
}
```

# Slides
* https://github.com/nakatamaho/mplapack/blob/v2.0/doc/presentation/20211128_%E7%B2%BE%E5%BA%A6%E4%BF%9D%E8%A8%BCmeeting.pdf

# How to build on Linux and Win (using Docker; recommended)

Ubuntu 20.04 (amd64, AArch64)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ git checkout v1.0.1
$ /usr/bin/time docker build -t mplapack:ubuntu2004 -f Dockerfile_ubuntu20.04 . 2>&1 | tee log.ubuntu2004
```

CentOS 7 (amd64, AArch64)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ git checkout v1.0.1
$ /usr/bin/time docker build -t mplapack:centos7 -f Dockerfile_CentOS7 . 2>&1 | tee log.CentOS7
```

CentOS 8 (amd64, AArch64)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ git checkout v1.0.1
$ /usr/bin/time docker build -t mplapack:centos8 -f Dockerfile_CentOS8 . 2>&1 | tee log.CentOS8
```

Ubuntu 18.04 (amd64, AArch64)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ git checkout v1.0.1
$ /usr/bin/time docker build -t mplapack:ubuntu1804 -f Dockerfile_ubuntu18.04 . 2>&1 | tee log.ubuntu1804
```

Ubuntu 20.04 (using Intel oneAPI)
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ git checkout v1.0.1
$ /usr/bin/time docker build -t mplapack:ubuntu2004intel -f Dockerfile_ubuntu20.04_intel . 2>&1 | tee log.ubuntu2004.intel
```

Windows 64bit (using cross compiler on Ubuntu) 
```
$ git clone https://github.com/nakatamaho/mplapack/
$ cd mplapack
$ git checkout v1.0.1
$ /usr/bin/time docker build -t mplapack:mingw -f  Dockerfile_ubuntu20.04_mingw64 . 2>&1 | tee log.mingw
```

# How to build and install (on Mac, using macports)

```
$ sudo port install gcc9 coreutils git ccache
$ git clone https://github.com/nakatamaho/mplapack.git --depth 1
$ cd mplapack
$ git checkout v1.0.1
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

# MPLAPACK Release schedule

## MPLAPACK 2.0.0 Release schedule
This is the release schedule for MPLAPACK 2.0.0
| Action | Date | Status | Description |
| --- | --- | --- | --- |
| complex lin (lapack routine; mpfr)      | 2021/11/29 | ng |Failed: CHR, CHK, CSR, CSK, CQ3, CTZ, CTS| 
| complex lin (driver routine; mpfr)      | 2021/11/29 | ng |Failed: CGB, CLS          | 
| complex lin (lapack routine; gmp)       | 2021/11/29 | ng |Failed: CSR, CSK, CQ3, CTZ|
| complex lin (driver routine; gmp)       | 2021/11/29 | ng |Falied: CGB               | 
| complex lin (lapack routine; dd)        | 2021/11/30 | ng |Failed: CGE CGB CGT CPO CPP CPB CPT CHE CHR CHK CH2 CSA CS2 CHP CSY CSR CSK CSP CTR CTP CTB CQR CRQ CLQ CQL CQ3 CTZ |
| complex lin (driver routine; dd)        | 2021/11/30 | ng |Failed: CGE CGB CGT CPO CPP CPB CHE CSA CH2 CHP CSY CSR CSK CSP CLS |
| complex lin (lapack routine; qd)        | 2021/11/30 | ng |Failed: CGE CGB CGT CPO CPP CPB CHE CHR CHK CSA CS2 CHP CSY CSR CSK CSP CTR CTP CTB CQR CRQ CLQ CQL CQ3 CTZ|
| complex lin (driver routine; qd)        | 2021/11/30 | ng |Failed: CGE CGB CGT CPO CPP  CPB  CHE  CSA CH2 CHP CSY CSR CSK CSP CLS |
| complex lin (lapack routine; _Float128) | 2021/11/29 | ng |Failed: CGE CGB CGT CPO CPP CPB CHE CHR CHK CHA CSA CHP CSY CSR CSK CSP CTR CTP CTB CQ3 CTZ|
| complex lin (driver routine; _Float128) | 2021/11/29 | ng |Failed:CGE CGB CGT CPO CPP CPB CHE CHA CSA CHP CSY CSR CSK CSP CLS|
| complex lin (lapack routine; _Float64x) | 2021/11/30 | ng |Failed: CHR CHK CSR CSK CQ3 CTZ CTS |
| complex lin (driver routine; _Float64x) | 2021/11/30 | ng |Failed: CGB CLS |
| complex lin (lapack routine; double)    | 2021/11/29 | ng |Failed: CHR CHK CSR CSK CQ3 CTZ CTQ CTS |
| complex lin (driver routine; double)    | 2021/11/29 | ng |Failed: CGB CLS         |
| sqrt fixes for QD and DD                | 2021/12/12 |     | near overflow |
| Impliment complex eig                |            |     |                          | 
| Impliment RFP version                |            |     |                          | 
| Impliment mixed precision version    |            |     |                          | 
| make tar ball for distribution       |            |     |                          |
| Add more examples                    |            |     |                          | 

## MPLAPACK 3.0.0 Release schedule
This is the release schedule for MPLAPACK 3.0.0
| Action | Date | Status | Description |
| --- | --- | --- | --- |
| Performance tuning using FAM on amd64          |            |     |                     |
| Fix Rlamch for QD, DD                          |            |     | svd.in fails.used with caution | 
| cleanup pow (REAL, long int)                   |            |     |                          | 
| Get rid of compiler warnings                   |            |     |                          | 
| lp64 ilp64 llp64 ilp32 lp32 cleanup            |            |     |                          | 
| FMA for QD, DD                                 |            |     |                     |
| Drop GMP version                               |            |     | Since trigonometric functions req'ed | 
| more benchmarks                                |            |     |                          | 
| add gmpfrxx                                    |            |     | https://math.berkeley.edu/~wilken/code/gmpfrxx/ | 
| optimized implementations                      |            |     |                          | 

## old schedules
* version 1.0.0 https://github.com/nakatamaho/mplapack/blob/master/doc/Release1.0.0.md

# History
* 2021/11/1 1.0.1 release. Fixing dd and qd arithmetic by Inte One API.
* 2021/10/1 1.0.0 release. Huge improvement; all real LAPACK routines are available; SVD, eigenproblem solver for non symmetric matrices are added. manual is available:  https://raw.githubusercontent.com/nakatamaho/mplapack/master/doc/manual/manual.pdf 
* 2021/4/11 0.9.3 release. CentOS7 AArch64 support
* 2021/4/6  0.9.1 release. CentOS support
* 2021/4/1  0.9.0 release. Rename to mplapack. You must rename include files, etc. Rewrite and recompilation required.
* 2012/12/25: MPACK 0.8.0 NVIDIA C2050 support for Rgemm in double-double, and preliminary support for Intel Xeon Phi. ANNOUNCE
* 2012/12/20: MPACK 0.8.0-RC2 Build fixes on various platforms.
* 2012/12/05: Our Rgemm dd paper "A Fast implementation of matrix-matrix product in double-double precision on NVIDIA C2050 and application to semidefinite programming" is selected as the Best Papers of The Third International Conference on Networking and Computing. Slide is here.
* 2012/11/29: MPACK 0.8.0-RC1 CUDA version of Rgemm in double-double precision is integrated.
* 2012/10/13: CUDA 4.2 or later version of accelerated Rgemm in double-double precision on NVIDIA C2050 GPU is now available. Note it does not work on CUDA 4.1. Original release announce is here., and preprint is available from here, and it will be presented at The Third International Conference on Networking and Computing Okinawa, Japan, December 5-7, 2012 .
* 2012/06/16: MPACK 0.7.0! Announcement
* 2012/06/16: Development has been migrated to SVN repository.
* 2011/10/28: Rgemm acceleration in double-double precision on NVIDIA C2050 GPU is now available. Even though papers are not published, you can just try by "make". Note that only CUDA 3.2 is supported. Origial release announce is here.
* 2011/08/24: Rgemm acceleration on NVIDIA C2050 GPU is coming. Unfortunately paper are rejected, so please wait... Here is a pdf slide.
* 2010/08/20: MPACK 0.6.7! Includes condition number estimators; Rgecon and Rpocon. Now 100 MLAPACK routines, and license has been changed to 2-clause BSD style license. No longer use LGPLv3.
* 2010/08/6: MPACK 0.6.6! Build fix release. Tested on various Linux distributions.
* 2010/05/31: A paper for MPACK (0.6.4) in Japanese has been uploaded.
* 2010/05/21: MPACK 0.6.5! MPFR support, and MBLAS license has been changed to BSD style. Still MLAPACK part is LGPLv3. I'll replace hopefully soon.
* 2010/01/13: MPACK 0.6.4! BUILD FIX RELEASE! PLEASE CHECK ON YOUR ENVIRONMENT! THANKS! ALSO WINDOWS IS SUPPORTED!
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

# Acknowledgment:

This work has been supported by:
The Special Postdoctoral Researchers' Program of RIKEN (2008, 2009)
Grant-in-Aid for Scientific Research (B) 21300017 from the Japan Society for the Promotion of Science (2009, 2010, 2011).
Microsoft Research CORE6 (2010), and the Japan Society for the Promotion of Science (JSPS KAKENHI Grant no. 18H03206) and TIS inc.

Also the M.N would like to thank Dr. Imamura Toshiyuki. Dr. Nakasato Naohito, Dr. Fujisawa Katsuki, Dr. Kouya Tomonori, Dr. Takahashi Daisuke, Dr. Goto Kazushige, Dr. Himeno Ryutaro, Dr. Hishimuna Toshiaki, Dr. Katagiri Takahiro, Dr. Ogita Takeshi, Dr. Kashiwagi Masahide, Dr. Yuasa Fukuko, Dr. Ishikawa Tadashi, Dr. Geshi Masaaki and Mr. Minato Yuichiro for warm encouragement.

# Citation
```
@misc{2109.13406,
Author = {Maho Nakata},
Title = {MPLAPACK version 1.0.0 user manual},
Year = {2021},
Eprint = {arXiv:2109.13406},
}
```

# contact
NAKATA Maho <maho.nakata@gmail.com> <maho@riken.jp>

