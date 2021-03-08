/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clatbs.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/

#include <mpblas.h>
#include <mplapack.h>

void
Clatbs(const char *uplo, const char *trans, const char *diag,
       const char *normin, INTEGER n, INTEGER kd, COMPLEX * AB, INTEGER ldab, COMPLEX * x, REAL * scale, REAL * cnorm, INTEGER * info)
{
    INTEGER i, j;
    REAL xj, rec, tjj;
    INTEGER jinc, jlen;
    REAL xbnd;
    INTEGER imax;
    REAL tmax;
    COMPLEX tjjs;
    REAL xmax, grow;
    INTEGER maind;
    REAL tscal;
    COMPLEX uscal;
    INTEGER jlast;
    COMPLEX csumj;
    INTEGER upper;
    REAL bignum;
    INTEGER notran;
    INTEGER jfirst;
    REAL smlnum;
    INTEGER nounit;
    REAL Half = .5, One = 1.0, Zero = 0.0, Two = 2.0;
    REAL mtemp1, mtemp2;

    *info = 0;
    upper = Mlsame(uplo, "U");
    notran = Mlsame(trans, "N");
    nounit = Mlsame(diag, "N");
//Test the input parameters.
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (!notran && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
	*info = -2;
    } else if (!nounit && !Mlsame(diag, "U")) {
	*info = -3;
    } else if (!Mlsame(normin, "Y") && !Mlsame(normin, "N")) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (kd < 0) {
	*info = -6;
    } else if (ldab < kd + 1) {
	*info = -8;
    }
    if (*info != 0) {
	Mxerbla("Clatbs", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Determine machine dependent parameters to control overflow.
    smlnum = Rlamch("Safe minimum");
    bignum = One / smlnum;
    smlnum = smlnum / Rlamch("Precision");
    bignum = One / smlnum;
    *scale = One;
    if (Mlsame(normin, "N")) {
//Compute the 1-norm of each column, not including the diagonal.
	if (upper) {
//A is upper triangular.
	    for (j = 0; j < n; j++) {
		jlen = min(kd, j - 1);
		cnorm[j] = RCasum(jlen, &AB[kd + 1 - jlen + j * ldab], 1);
	    }
	} else {
//A is lower triangular.
	    for (j = 0; j < n; j++) {
		jlen = min(kd, n - j);
		if (jlen > 0) {
		    cnorm[j] = RCasum(jlen, &AB[j * ldab + 2], 1);
		} else {
		    cnorm[j] = Zero;
		}
	    }
	}
    }
//Scale the column norms by TSCAL if the maximum element in CNORM is
//greater than BIGNUM/Two
    imax = iRamax(n, &cnorm[1], 1);
    tmax = cnorm[imax];
    if (tmax <= bignum * Half) {
	tscal = Zero;
    } else {
	tscal = Half / (smlnum * tmax);
	Rscal(n, tscal, &cnorm[1], 1);
    }
//Compute a bound on the computed solution vector to see if the
//Level 2 BLAS routine ZTBSV can be used.
    xmax = Zero;
    for (j = 0; j < n; j++) {
	mtemp1 = xmax, mtemp2 = abs(x[j].real() / Two) + abs(x[j].imag() / Two);
	xmax = max(mtemp1, mtemp2);
    }
    xbnd = xmax;
    if (notran) {
//Compute the growth in A * x = b.
	if (upper) {
	    jfirst = n;
	    jlast = 1;
	    jinc = -1;
	    maind = kd + 1;
	} else {
	    jfirst = 1;
	    jlast = n;
	    jinc = 1;
	    maind = 1;
	}

	if (tscal != One) {
	    grow = Zero;
	    goto L60;
	}

	if (nounit) {
//A is non-unit triangular.
//Compute GROW = 1/G(j) and XBND = 1/M(j).
//Initially, G(0) = max{x(i), i=1,...,n}.
	    grow = Half / max(xbnd, smlnum);
	    xbnd = grow;
	    for (j = jfirst; j <= jlast; j = j + jinc) {
//Exit the loop if the growth factor is too small.
		if (grow <= smlnum) {
		    goto L60;
		}
		tjjs = AB[maind + j * ldab];
		tjj = abs(tjjs.real()) + abs(tjjs.imag());
		if (tjj >= smlnum) {
//M(j) = G(j-1) / abs(A(j,j))
		    mtemp1 = xbnd, mtemp2 = min(One, tjj) * grow;
		    xbnd = min(mtemp1, mtemp2);
		} else {
//M(j) could overflow, set XBND to Zero
		    xbnd = Zero;
		}
		if (tjj + cnorm[j] >= smlnum) {
//G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
		    grow *= tjj / (tjj + cnorm[j]);
		} else {
//G(j) could overflow, set GROW to Zero
		    grow = Zero;
		}
	    }
	    grow = xbnd;
	} else {
//A is unit triangular.
//Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
	    mtemp1 = One, mtemp2 = Half / max(xbnd, smlnum);
	    grow = min(mtemp1, mtemp2);
	    for (j = jfirst; j <= jlast; j = j + jinc) {
//Exit the loop if the growth factor is too small.
		if (grow <= smlnum) {
		    goto L60;
		}
//G(j) = G(j-1)*( 1 + CNORM(j) )
		grow = grow * One / (cnorm[j] + One);
	    }
	}
      L60:
	;
    } else {
//Compute the growth in A**T * x = b  or  A**H * x = b.
	if (upper) {
	    jfirst = 1;
	    jlast = n;
	    jinc = 1;
	    maind = kd + 1;
	} else {
	    jfirst = n;
	    jlast = 1;
	    jinc = -1;
	    maind = 1;
	}

	if (tscal != One) {
	    grow = Zero;
	    goto L90;
	}

	if (nounit) {
//A is non-unit triangular.
//Compute GROW = 1/G(j) and XBND = 1/M(j).
//Initially, M(0) = max{x(i), i=1,...,n}.
	    grow = Half / max(xbnd, smlnum);
	    xbnd = grow;
	    for (j = jfirst; j <= jlast; j = j + jinc) {
//Exit the loop if the growth factor is too small.
		if (grow <= smlnum) {
		    goto L90;
		}
//G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
		xj = cnorm[j] + One;
		mtemp1 = grow, mtemp2 = xbnd / xj;
		grow = min(mtemp1, mtemp2);
		tjjs = AB[maind + j * ldab];
		tjj = abs(tjjs.real()) + abs(tjjs.imag());
		if (tjj >= smlnum) {
//M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
		    if (xj > tjj) {
			xbnd = xbnd * tjj / xj;
		    }
		} else {
//M(j) could overflow, set XBND to Zero
		    xbnd = Zero;
		}
	    }
	    grow = min(grow, xbnd);
	} else {
//A is unit triangular.
//Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
	    mtemp1 = One, mtemp2 = Half / max(xbnd, smlnum);
	    grow = min(mtemp1, mtemp2);
	    for (j = jfirst; j <= jlast; j += jinc) {
//Exit the loop if the growth factor is too small.
		if (grow <= smlnum) {
		    goto L90;
		}
// G(j) = ( 1 + CNORM(j) )*G(j-1)
		xj = cnorm[j] + One;
		grow = grow / xj;
	    }
	}
      L90:
	;
    }
    if (grow * tscal > smlnum) {
//Use the Level 2 BLAS solve if the reciprocal of the bound on
//elements of X is not too small.
	Ctbsv(uplo, trans, diag, n, kd, &AB[0], ldab, &x[0], 1);
    } else {
//Use a Level 1 BLAS solve, scaling intermediate results.
	if (xmax > bignum * Half) {
//Scale X so that its components are less than or equal to
//BIGNUM in absolute value.
	    *scale = bignum * Half / xmax;
	    CRscal(n, *scale, &x[0], 1);
	    xmax = bignum;
	} else {
	    xmax = xmax * Two;
	}

	if (notran) {
//Solve A * x = b
	    for (j = jfirst; j <= jlast; j = j + jinc) {
//Compute x(j) = b(j) / A(j,j), scaling x if necessary.
		xj = abs(x[j].real()) + abs(x[j].imag());
		if (nounit) {
		    tjjs = AB[maind + j * ldab] * tscal;
		} else {
		    tjjs = tscal;
		    if (tscal == One) {
			goto L110;
		    }
		}
		tjj = abs(tjjs.real()) + abs(tjjs.imag());
		if (tjj > smlnum) {
//abs(A(j,j)) > SMLNUM:
		    if (tjj < One) {
			if (xj > tjj * bignum) {
//Scale x by 1/b(j).
			    rec = One / xj;
			    CRscal(n, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
		    }
		    x[j] = Cladiv(x[j], tjjs);
		    xj = abs(x[j].real()) + abs(x[j].imag());
		} else if (tjj > Zero) {
//0 < abs(A(j,j)) <= SMLNUM:
		    if (xj > tjj * bignum) {
//Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
//to avoid overflow when dividing by A(j,j).
			rec = tjj * bignum / xj;
			if (cnorm[j] > One) {
//Scale by 1/CNORM(j) to avoid overflow when
//multiplying x(j) times column j.
			    rec = rec / cnorm[j];
			}
			CRscal(n, rec, &x[0], 1);
			*scale = *scale * rec;
			xmax = xmax * rec;
		    }
		    x[j] = Cladiv(x[j], tjjs);
		    xj = abs(x[j].real()) + abs(x[j].imag());
		} else {
//A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
//scale = 0, and compute a solution to A*x = Zero
		    for (i = 0; i < n; i++) {
			x[i] = Zero;
		    }
		    x[j] = One;
		    xj = Zero;
		    *scale = Zero;
		    xmax = Zero;
		}
	      L110:
//Scale x if necessary to avoid overflow when adding a
//multiple of column j of A.
		if (xj > One) {
		    rec = One / xj;
		    if (cnorm[j] > (bignum - xmax) * rec) {
//Scale x by 1/(2*abs(x(j))).
			rec = rec * Half;
			CRscal(n, rec, &x[0], 1);
			*scale = *scale * rec;
		    }
		} else if (xj * cnorm[j] > bignum - xmax) {
//Scale x by 1/Two
		    CRscal(n, Two, &x[0], 1);
		    *scale *= Half;
		}

		if (upper) {
		    if (j > 1) {
//Compute the update
//   x(max((INTEGER)1,j-kd):j-1) := x(max((INTEGER)1,j-kd):j-1) -
//                         x(j)* A(max((INTEGER)1,j-kd):j-1,j)
			jlen = min(kd, j - 1);
			Caxpy(jlen, -x[j] * tscal, &AB[kd + 1 - jlen + j * ldab], 1, &x[j - jlen], 1);
			i = iCamax(j - 1, &x[0], 1);
			xmax = abs(x[i].real()) + abs(x[i].imag());
		    }
		} else if (j < n) {
//Compute the update
//   x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
//                         x(j) * A(j+1:min(j+kd,n),j)
		    jlen = min(kd, n - j);
		    if (jlen > 0) {
			Caxpy(jlen, -x[j] * tscal, &AB[j * ldab + 2], 1, &x[j + 1], 1);
		    }
		    i = j + iCamax(n - j, &x[j + 1], 1);
		    xmax = abs(x[i].real()) + abs(x[i].imag());
		}
	    }
	} else if (Mlsame(trans, "T")) {
//Solve A**T * x = b
	    for (j = jfirst; j <= jlast; j = j + jinc) {
//Compute x(j) = b(j) - sum A(k,j)*x(k).
//                      k<>j
		xj = abs(x[j].real()) + abs(x[j].imag());
		uscal = tscal;
		rec = One / max(xmax, One);
		if (cnorm[j] > (bignum - xj) * rec) {
//If x(j) could overflow, scale x by 1/(2*XMAX).
		    rec = rec * Half;
		    if (nounit) {
			tjjs = AB[maind + j * ldab] * tscal;
		    } else {
			tjjs = tscal;
		    }
		    tjj = abs(tjjs.real()) + abs(tjjs.imag());
		    if (tjj > One) {
//Divide by A(j,j) when scaling x if A(j,j) > One
			mtemp1 = One, mtemp2 = rec * tjj;
			rec = min(mtemp1, mtemp2);
			uscal = Cladiv(uscal, tjjs);
		    }
		    if (rec < One) {
			CRscal(n, rec, &x[0], 1);
			*scale = *scale * rec;
			xmax = xmax * rec;
		    }
		}
		csumj = Zero;
		if (uscal == One) {
//If the scaling needed for A in the dot product is 1,
//call ZDOTU to perform the dot product.
		    if (upper) {
			jlen = min(kd, j - 1);
			csumj = Cdotu(jlen, &AB[kd + 1 - jlen + j * ldab], 1, &x[j - jlen], 1);
		    } else {
			jlen = min(kd, n - j);
			if (jlen > 1) {
			    csumj = Cdotu(jlen, &AB[j * ldab + 2], 1, &x[j + 1], 1);
			}
		    }
		} else {
//Otherwise, use in-line code for the dot product.
		    if (upper) {
			jlen = min(kd, j - 1);
			for (i = 0; i < jlen; i++) {
			    csumj = csumj + AB[kd + i - jlen + j * ldab] * uscal * x[j - jlen - 1 + i];
			}
		    } else {
			jlen = min(kd, n - j);
			for (i = 0; i < jlen; i++) {
			    csumj = csumj + AB[i + 1 + j * ldab] * uscal * x[j + i];
			}
		    }
		}
		if (uscal.real() == tscal && uscal.imag() == Zero) {
//Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
//was not used to scale the dotproduct.
		    x[j] = x[j] - csumj;
		    xj = abs(x[j].real()) + abs(x[j].imag());
		    if (nounit) {
//Compute x(j) = x(j) / A(j,j), scaling if necessary.
			tjjs = AB[maind + j * ldab] * tscal;
		    } else {
			tjjs = tscal;
			if (tscal == One) {
			    goto L160;
			}
		    }
		    tjj = abs(tjjs.real()) + abs(tjjs.imag());
		    if (tjj > smlnum) {
//abs(A(j,j)) > SMLNUM:
			if (tjj < One) {
			    if (xj > tjj * bignum) {
//Scale X by 1/abs(x(j)).
				rec = One / xj;
				CRscal(n, rec, &x[0], 1);
				*scale = *scale * rec;
				xmax = xmax * rec;
			    }
			}
			x[j] = Cladiv(x[j], tjjs);
		    } else if (tjj > Zero) {
//0 < abs(A(j,j)) <= SMLNUM:
			if (xj > tjj * bignum) {
//Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
			    rec = tjj * bignum / xj;
			    CRscal(n, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
			x[j] = Cladiv(x[j], tjjs);
		    } else {
//A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
//scale = 0 and compute a solution to A**T *x = Zero
			for (i = 0; i < n; i++) {
			    x[i] = Zero;
			}
			x[j] = One;
			*scale = Zero;
			xmax = Zero;
		    }
		  L160:
		    ;
		} else {
//Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
//product has already been divided by 1/A(j,j).
		    x[j] = Cladiv(x[j], tjjs) - csumj;
		}
		mtemp1 = xmax, mtemp2 = abs(x[j].real()) + abs(x[j].imag());
		xmax = max(mtemp1, mtemp2);
	    }
	} else {
//Solve A**H * x = b
	    for (j = jfirst; j <= jlast; j += jinc) {
//Compute x(j) = b(j) - sum A(k,j)*x(k).
//                      k<>j
		xj = abs(x[j].real()) + abs(x[j].imag());
		uscal = tscal;
		rec = One / max(xmax, One);
		if (cnorm[j] > (bignum - xj) * rec) {
//If x(j) could overflow, scale x by 1/(2*XMAX).
		    rec = rec * Half;
		    if (nounit) {
			tjjs = conj(AB[maind + j * ldab]) * tscal;
		    } else {
			tjjs = tscal;
		    }
		    tjj = abs(tjjs.real()) + abs(tjjs.imag());
		    if (tjj > One) {
//Divide by A(j,j) when scaling x if A(j,j) > One
			mtemp1 = One, mtemp2 = rec * tjj;
			rec = min(mtemp1, mtemp2);
			uscal = Cladiv(uscal, tjjs);
		    }
		    if (rec < One) {
			CRscal(n, rec, &x[0], 1);
			*scale = *scale * rec;
			xmax = xmax * rec;
		    }
		}
		csumj = Zero;
		if (uscal == One) {
//If the scaling needed for A in the dot product is 1,
//call ZDOTC to perform the dot product.
		    if (upper) {
			jlen = min(kd, j - 1);
			csumj = Cdotc(jlen, &AB[kd + 1 - jlen + j * ldab], 1, &x[j - jlen], 1);
		    } else {
			jlen = min(kd, n - j);
			if (jlen > 1) {
			    csumj = Cdotc(jlen, &AB[j * ldab + 2], 1, &x[j + 1], 1);
			}
		    }
		} else {
//Otherwise, use in-line code for the dot product.
		    if (upper) {
			jlen = min(kd, j - 1);
			for (i = 0; i < jlen; i++) {
			    csumj = csumj + conj(AB[kd + i - jlen + j * ldab]) * uscal * x[j - jlen - 1 + i];
			}
		    } else {
			jlen = min(kd, n - j);
			for (i = 0; i < jlen; i++) {
			    csumj = csumj + conj(AB[i + 1 + j * ldab]) * uscal * x[j + i];
			}
		    }
		}
		if (uscal.real() == tscal && uscal.imag() == Zero) {
//Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
//was not used to scale the dotproduct.
		    x[j] = x[j] - csumj;
		    xj = abs(x[j].real()) + abs(x[j].imag());
		    if (nounit) {
//Compute x(j) = x(j) / A(j,j), scaling if necessary.
			tjjs = conj(AB[maind + j * ldab]) * tscal;
		    } else {
			tjjs = tscal;
			if (tscal == One) {
			    goto L210;
			}
		    }
		    tjj = abs(tjjs.real()) + abs(tjjs.imag());
		    if (tjj > smlnum) {
//abs(A(j,j)) > SMLNUM:
			if (tjj < One) {
			    if (xj > tjj * bignum) {
//Scale X by 1/abs(x(j)).
				rec = One / xj;
				CRscal(n, rec, &x[0], 1);
				*scale = *scale * rec;
				xmax = xmax * rec;
			    }
			}
			x[j] = Cladiv(x[j], tjjs);
		    } else if (tjj > Zero) {
//0 < abs(A(j,j)) <= SMLNUM:
			if (xj > tjj * bignum) {
//Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
			    rec = tjj * bignum / xj;
			    CRscal(n, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
			x[j] = Cladiv(x[j], tjjs);
		    } else {
//A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
//scale = 0 and compute a solution to A**H *x = Zero
			for (i = 0; i < n; i++) {
			    x[i] = Zero;
			}
			x[j] = One;
			*scale = Zero;
			xmax = Zero;
		    }
		  L210:
		    ;
		} else {
//Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
//product has already been divided by 1/A(j,j).
		    x[j] = Cladiv(x[j], tjjs) - csumj;
		}
		mtemp1 = xmax, mtemp2 = abs(x[j].real()) + abs(x[j].imag());
		xmax = max(mtemp1, mtemp2);
	    }
	}
	*scale = *scale / tscal;
    }
//Scale the column norms by 1/TSCAL for return.
    if (tscal != One) {
	Rscal(n, One / tscal, &cnorm[1], 1);
    }
    return;
}
