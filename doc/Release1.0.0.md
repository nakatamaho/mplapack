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
| Build on Ubuntu 20.04 amd64          | 2021-09-27 |1e88da2| _Float128: supp. by libc only binary128, _Float64x == long double |
| Build on Ubuntu 18.04 amd64          | 2021-09-27 |1e88da2| _Float128: supp. by libc only binary128, _Float64x == long double |
| Build on Ubuntu 20.04 Intel oneAPI   | 2021-09-27 |1e88da2| _Float128: supp. by libc only binary128, _Float64x == long double |
| Build on Ubuntu 20.04 mingw64        | 2021-09-27 |1e88da2| _Float128: supp. by libquadmath, _Float64x == __float80 == long double |
| Build on CentOS7 amd64               | 2021-09-27 |1e88da2| _Float128: supp. by libquadmath, _Float64x == __float80 == long double |
| Build on CentOS8 amd64               | 2021-09-27 |1e88da2| _Float128: supp. by libc only binary128, _Float64x == long double |
| Build on MacOS Big Sur amd64         | 2021-09-27 |1e88da2| _Float128: supp. by libquadmath, _Float64x == __float80 == long double |
| Build on Ubuntu 20.04 AArch64        | 2021-09-27 |1e88da2| _Float128: supp. by libc and _Float128 == long double |
| Build on CentOS7 AArch64             | 2021-09-27 |1e88da2| long double == binary128, _Float128 nor __float128 are supported |
| Build on CentOS8 AArch64             | 2021-09-27 |1e88da2| _Float128: supp. by libc and _Float128 == long double |
| Build on Ubuntu 20.04 PPC64LE        | 2021-09-27 |1e88da2| _Float128: supp. by libc and _Float128 != long double |
| Add documents                        | 2021-09-25 |       |                          | 
| Upload to arXiv                      | 2021-09-27 |       | https://arxiv.org/abs/2109.13406 | 
| Add bibtex citation                  | 2021-09-29 |       |                          | 
| add tag                              | 2021-09-27 |       |                          | 
| move tag after appearing at arXiv    | 2021-09-29 |       |                          | 
| Release                              | 2021-10-01 |       |                          | 

