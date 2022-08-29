## MPLAPACK 2.0.0 Release process
This is the release process for MPLAPACK 2.0.0
| Action | Date | Status | Description |
| --- | --- | --- | --- |
| complex lin (lapack routine; mpfr)      | 2022/5/6, 2022/03/21 | ok |               | 
| complex lin (driver routine; mpfr)      | 2022/5/6, 2022/03/21 | ok |               | 
| complex lin (lapack routine; gmp)       | 2022/5/6, 2022/03/21 | ok |               |
| complex lin (driver routine; gmp)       | 2022/5/6, 2022/03/21 | ok |               | 
| complex lin (lapack routine; dd)        | 2022/5/6, 2022/03/21 | ok |               | 
| complex lin (driver routine; dd)        | 2022/5/6, 2022/03/21 | ok |               |
| complex lin (lapack routine; qd)        | 2022/05/09, 2022/03/21 | ok |        | 
| complex lin (driver routine; qd)        | 2022/05/09, 2022/03/21 |ok | |
| complex lin (lapack routine; _Float128) | 2022/5/6, 2022/03/21 | ok |  |
| complex lin (driver routine; _Float128) | 2022/5/6, 2022/03/21 | ok |  |
| complex lin (lapack routine; _Float64x) | 2022/5/6, 2022/03/21 | ok |  |
| complex lin (driver routine; _Float64x) | 2022/5/6, 2022/03/21 | ok |  |
| complex lin (lapack routine; double)    | 2022/5/6, 2022/03/21 | ok |  |
| complex lin (driver routine; double)    | 2022/5/6, 2022/03/21 | ok |  |
| sqrt fixes for QD and DD                | 2021/12/12 |     | near overflow |
| complex eig Cbal.in                     | 2022/03/25 | ok  | Testing the balancing of a general matrix  | 
| complex eig Cbak.in                     | 2022/03/26 | ok  | Testing backward transformation | 
| complex eig Cgbal.in                    | 2022/03/27 | ok  | Testing the balancing of a general matrix  | 
| complex eig Cgbak.in                    | 2022/04/21 | ok  | Testing the back transformation of a pair of COMPLEX balanced matrices | 
| complex eig Cec.in                      | 2022/04/24 | ok  | Testing COMPLEX Eigen Condition Routines |
| complex eig nep.in                      | 2022/04/26 | ok  | testing Nonsymmetric Eigenvalue Problem routines  for COMPLEX |
| complex eig lse.in                      | 2022/04/26 | ok  | testing Constrained Linear Least Squares routines for COMPLEX |
| complex eig Csb.in                      | 2022/04/26 | ok  | testing Hermitian Eigenvalue Problem routines                 |
| complex eig Cglm.in                     | 2022/04/26 | ok  | Tests of the Generalized Linear Regression Model routines     |
| complex eig Cgqr.in                     | 2022/04/26 | ok  | testing Generalized QR and RQ routines |
| complex eig Cbb.in                      | 2022/04/26 | ok  | reduction of a general band matrix to real bidiagonal form |
| complex eig sep.in                      | 2022/04/30 | ok  | testing Hermitian Eigenvalue Problem routines for complex |
| complex eig se2.in                      | 2022/05/01 | ok  | testing Hermitian Eigenvalue Problem routines for complex (two stages ver.) |
| complex eig Ced.in                      | 2022/05/02 | ok  | testing Hermitian Eigenvalue Problem Expert drivers and drivers for complex |
| complex eig csd.in                      | 2022/05/10 | ok  | testing complex CS decomposition routines       |
| complex eig Cgg.in                      | 2022/05/04 | ok | testing Nonsymmetric Eigenvalue Problem routines  |
| complex eig Cgd.in                      | 2022/05/08 | ok | testing Complex Nonsymmetric generalized Eigenvalue/Schur Form Driver/Expert Driver |
| complex eig Csg.in                      | 2022/05/09 | ok | Tests of the complex Generalized Hermitian Eigenvalue Problem routines |
| complex eig svd.in                      | 2022/05/06 | NG for dd, qd | Testing Singular Value Decomposition routines for complex matrices. In some cases, Cgesvj has large error for dd and qd. Othre than that, it is ok|
| complex eig gsv.in                      | 2022/05/09 |  ok   | testing complex Generalized SVD routines |
| Impliment RFP version                   | 2022/05/16 |  ok   |  Test for linear equation routines RFP format (real, complex) | 
| libqd overflow fix                      | 2022/07/23 |  ok   |  sqrt(Rlamch("O") and Rlamch("O")/1.0 failed | 
| Fixes for Rgejsv, Cgejsv, Rgesvj and Cgesvj using dd/qd | 2022/07/24 | ok |  |
| Add more examples                       | 2022/06/14 |     |                          | 
| make tar ball for distribution          | 2022/06/14 |     |                          |
| alpha release                           | 2022/06/14 |     |                          |
| update document                         |            |     |                          |
| release                                 | 2022/7/26  |     |                          |

