/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clarfb.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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
 * Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.
 * 
 * $COPYRIGHT$
 * 
 * Additional copyrights may follow
 * 
 * $HEADER$
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * - Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer listed in this
 * license in the documentation and/or other materials provided with the
 * distribution.
 * 
 * - Neither the name of the copyright holders nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <mblas.h>
#include <mlapack.h>

void
Clarfb(const char *side, const char *trans, const char *direct,
       const char *storev, INTEGER m, INTEGER n, INTEGER k, COMPLEX * V, INTEGER ldv, COMPLEX * T, INTEGER ldt, COMPLEX * C, INTEGER ldc, COMPLEX * work, INTEGER ldwork)
{
    INTEGER i, j;
    REAL One = 1.0;
    char transt;

    //Quick return if possible
    if (m <= 0 || n <= 0) {
	return;
    }
    if (Mlsame(trans, "N")) {
	transt = 'C';
    } else {
	transt = 'N';
    }
    if (Mlsame(storev, "C")) {
	if (Mlsame(direct, "F")) {
//Let V = (V1) (first K rows)
//        (V2)
//where V1 is unit lower triangular.
	    if (Mlsame(side, "L")) {
//Form H * C or H' * C  where  C = ( C1 )
//                                 ( C2 )
//W: = C' * V  =  (C1' * V1 + C2'*V2)  (stored in WORK)
//W: = C1'
		for (j = 0; j < k; j++) {
		    Ccopy(n, &C[j], ldc, &work[j * ldwork], 1);
		    Clacgv(n, &work[j * ldwork], 1);
		}
//W: = W * V1
		Ctrmm("Right", "Lower", "No transpose", "Unit", n, k, (COMPLEX) One, V, ldv, work, ldwork);
		if (m > k) {
//W: = W + C2'*V2
		    Cgemm("Conjugate transpose", "No transpose", n, k, m - k, (COMPLEX) One, &C[k], ldc, &V[k], ldv, (COMPLEX) One, work, ldwork);
		}
//W: = W * T'  or  W * T
		Ctrmm("Right", "Upper", &transt, "Non-unit", n, k, (COMPLEX) One, T, ldt, work, ldwork);
//C: = C - V * W'
		if (m > k) {
//C2: = C2 - V2 * W'
		    Cgemm("No transpose", "Conjugate transpose", m - k, n, k, (COMPLEX) - One, &V[k], ldv, work, ldwork, (COMPLEX) One, &C[k], ldc);
		}
//W: = W * V1'
		Ctrmm("Right", "Lower", "Conjugate transpose", "Unit", n, k, (COMPLEX) One, V, ldv, work, ldwork);
//C1: = C1 - W'
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[j + i * ldc] = C[j + i * ldc] - conj(work[i + j * ldwork]);
		    }
		}
	    } else if (Mlsame(side, "R")) {
//Form C * H or C * H'  where  C = ( C1  C2 )
//W: = C * V = (C1 * V1 + C2 * V2) (stored in WORK)
//W: = C1
		for (j = 0; j < k; j++) {
		    Ccopy(m, &C[j * ldc], 1, &work[j * ldwork], 1);
		}
//W: = W * V1
		Ctrmm("Right", "Lower", "No transpose", "Unit", m, k, (COMPLEX) One, V, ldv, work, ldwork);
		if (n > k) {
//W: = W + C2 * V2
		    Cgemm("No transpose", "No transpose", m, k, n - k, (COMPLEX) One, &C[k * ldc], ldc, &V[k], ldv, (COMPLEX) One, work, ldwork);
		}
//W: = W * T or W * T'
		Ctrmm("Right", "Upper", trans, "Non-unit", m, k, (COMPLEX) One, T, ldt, work, ldwork);
//C: = C - W * V'
		if (n > k) {
//C2: = C2 - W * V2'
		    Cgemm("No transpose", "Conjugate transpose", m, n - k, k, -(COMPLEX) One, work, ldwork, &V[k], ldv, (COMPLEX) One, &C[k * ldc], ldc);
		}
//W: = W * V1'
		Ctrmm("Right", "Lower", "Conjugate transpose", "Unit", m, k, (COMPLEX) One, V, ldv, work, ldwork);
//C1: = C1 - W
		for (j = 0; j < k; j++) {
		    for (i = 0; i < m; i++) {
			C[i + j * ldc] = C[i + j * ldc] - work[i + j * ldwork];
		    }
		}
	    }
	} else {
//Let V = (V1)
//        (V2) (last K rows)
//where V2 is unit upper triangular.
	    if (Mlsame(side, "L")) {
//Form H * C or H' * C  where  C = ( C1 )
//                                 ( C2 )
//W: = C' * V  =  (C1' * V1 + C2'*V2)  (stored in WORK)
//W: = C2'
		for (j = 0; j < k; j++) {
		    Ccopy(n, &C[m - k + j], ldc, &work[j * ldwork], 1);
		    Clacgv(n, &work[j * ldwork], 1);
		}
//W: = W * V2
		Ctrmm("Right", "Upper", "No transpose", "Unit", n, k, (COMPLEX) One, &V[m - k], ldv, work, ldwork);
		if (m > k) {
//W: = W + C1'*V1
		    Cgemm("Conjugate transpose", "No transpose", n, k, m - k, (COMPLEX) One, C, ldc, V, ldv, (COMPLEX) One, work, ldwork);
		}
//W: = W * T'  or  W * T
		Ctrmm("Right", "Lower", &transt, "Non-unit", n, k, (COMPLEX) One, T, ldt, work, ldwork);
//C: = C - V * W'
		if (m > k) {
//C1:= C1 - V1 * W'
		    Cgemm("No transpose", "Conjugate transpose", m - k, n, k, (COMPLEX) - One, V, ldv, work, ldwork, (COMPLEX) One, C, ldc);
		}
//W: = W * V2'
		Ctrmm("Right", "Upper", "Conjugate transpose", "Unit", n, k, (COMPLEX) One, &V[m - k], ldv, work, ldwork);
//C2:= C2 - W'
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[m - k + j + i * ldc] = C[m - k + j + i * ldc] - conj(work[i + j * ldwork]);
		    }
		}
	    } else if (Mlsame(side, "R")) {
//Form C * H or C * H'  where  C = ( C1  C2 )
//W: = C * V = (C1 * V1 + C2 * V2) (stored in WORK)
//W: = C2
		for (j = 0; j < k; j++) {
		    Ccopy(m, &C[(n - k + j) * ldc], 1, &work[j * ldwork], 1);
		}
//W:= W * V2
		Ctrmm("Right", "Upper", "No transpose", "Unit", m, k, (COMPLEX) One, &V[n - k], ldv, work, ldwork);
		if (n > k) {
//W:= W + C1 * V1
		    Cgemm("No transpose", "No transpose", m, k, n - k, (COMPLEX) One, C, ldc, V, ldv, (COMPLEX) One, work, ldwork);
		}
//W:= W * T or W * T
		Ctrmm("Right", "Lower", trans, "Non-unit", m, k, (COMPLEX) One, T, ldt, work, ldwork);
//C:= C - W * V'
		if (n > k) {
//C1:= C1 - W * V1'
		    Cgemm("No transpose", "Conjugate transpose", m, n - k, k, (COMPLEX) - One, work, ldwork, V, ldv, (COMPLEX) One, C, ldc);
		}
//W: = W * V2'
		Ctrmm("Right", "Upper", "Conjugate transpose", "Unit", m, k, (COMPLEX) One, &V[n - k], ldv, work, ldwork);
//C2:= C2 - W
		for (j = 0; j < k; j++) {
		    for (i = 0; i < m; i++) {
			C[i + (n - k + j) * ldc] = C[i + (n - k + j) * ldc] - work[i + j * ldwork];
		    }
		}
	    }
	}
    } else if (Mlsame(storev, "R")) {
	if (Mlsame(direct, "F")) {
//Let V = (V1 V2) (V1:first K columns)
//where V1 is unit upper triangular.

	    if (Mlsame(side, "L")) {
//Form H * C or H' * C  where  C = ( C1 )
//                                 ( C2 )
//W:= C' * V' = (C1'*V1' + C2'*V2') (stored in WORK)
//W:= C1'
		for (j = 0; j < k; j++) {
		    Ccopy(n, &C[j], ldc, &work[j * ldwork], 1);
		    Clacgv(n, &work[j * ldwork], 1);
		}
//W:= W * V1'
		Ctrmm("Right", "Upper", "Conjugate transpose", "Unit", n, k, (COMPLEX) One, V, ldv, work, ldwork);
		if (m > k) {
//W:= W + C2'*V2'
		    Cgemm("Conjugate transpose", "Conjugate transpose", n, k, m - k, (COMPLEX) One, &C[k], ldc, &V[k * ldv], ldv, (COMPLEX) One, work, ldwork);
		}
//W:= W * T'  or  W * T
		Ctrmm("Right", "Upper", &transt, "Non-unit", n, k, (COMPLEX) One, T, ldt, work, ldwork);
//C:= C - V' * W'
		if (m > k) {
//C2:= C2 - V2' * W'
		    Cgemm("Conjugate transpose", "Conjugate transpose", m - k, n, k, (COMPLEX) - One, &V[k * ldv], ldv, work, ldwork, (COMPLEX) One, &C[k], ldc);
		}
//W:= W * V1
		Ctrmm("Right", "Upper", "No transpose", "Unit", n, k, (COMPLEX) One, V, ldv, work, ldwork);
//C1:= C1 - W'
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[j + i * ldc] = C[j + i * ldc] - conj(work[i + j * ldwork]);
		    }
		}
	    } else if (Mlsame(side, "R")) {
//Form C * H or C * H'  where  C = ( C1  C2 )
//W:= C * V'  =  (C1*V1' + C2 * V2')  (stored in WORK)
//W:= C1
		for (j = 0; j < k; j++) {
		    Ccopy(m, &C[j * ldc], 1, &work[j * ldwork], 1);
		}
//W:= W * V1'
		Ctrmm("Right", "Upper", "Conjugate transpose", "Unit", m, k, (COMPLEX) One, V, ldv, work, ldwork);
		if (n > k) {
//W:= W + C2 * V2'
		    Cgemm("No transpose", "Conjugate transpose", m, k, n - k, (COMPLEX) One, &C[k * ldc], ldc, &V[k * ldv], ldv, (COMPLEX) One, work, ldwork);
		}
//W:= W * T or W * T'
		Ctrmm("Right", "Upper", trans, "Non-unit", m, k, (COMPLEX) One, T, ldt, work, ldwork);
//C:= C - W * V
		if (n > k) {
//C2:= C2 - W * V2
		    Cgemm("No transpose", "No transpose", m, n - k, k, -(COMPLEX) One, work, ldwork, &V[k * ldv], ldv, (COMPLEX) One, &C[k * ldc], ldc);
		}
//W:= W * V1
		Ctrmm("Right", "Upper", "No transpose", "Unit", m, k, (COMPLEX) One, V, ldv, work, ldwork);
//C1:= C1 - W
		for (j = 0; j < k; j++) {
		    for (i = 0; i < m; i++) {
			C[i + j * ldc] = C[i + j * ldc] - work[i + j * ldwork];
		    }
		}
	    }
	} else {
//Let V = (V1 V2) (V2:last K columns)
// where V2 is unit lower triangular.
	    if (Mlsame(side, "L")) {
//Form H * C or H' * C  where  C = ( C1 )
//                                 ( C2 )
//W:= C' * V' = (C1'*V1' + C2'*V2') (stored in WORK)
//W:= C2'
		for (j = 0; j < k; j++) {
		    Ccopy(n, &C[m - k + j], ldc, &work[j * ldwork], 1);
		    Clacgv(n, &work[j * ldwork], 1);
		}
//W:= W * V2'
		Ctrmm("Right", "Lower", "Conjugate transpose", "Unit", n, k, (COMPLEX) One, &V[(m - k) * ldv], ldv, work, ldwork);
		if (m > k) {
//W:= W + C1'*V1'
		    Cgemm("Conjugate transpose", "Conjugate transpose", n, k, m - k, (COMPLEX) One, C, ldc, V, ldv, (COMPLEX) One, work, ldwork);
		}
//W:= W * T'  or  W * T
		Ctrmm("Right", "Lower", &transt, "Non-unit", n, k, (COMPLEX) One, T, ldt, work, ldwork);
//C:= C - V' * W'
		if (m > k) {
//C1:= C1 - V1' * W'
		    Cgemm("Conjugate transpose", "Conjugate transpose", m - k, n, k, (COMPLEX) - One, V, ldv, work, ldwork, (COMPLEX) One, C, ldc);
		}
//W:= W * V2
		Ctrmm("Right", "Lower", "No transpose", "Unit", n, k, (COMPLEX) One, &V[(m - k) * ldv], ldv, work, ldwork);
//C2:= C2 - W'
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[m - k + j + i * ldc] = C[m - k + j + i * ldc] - conj(work[i + j * ldwork]);
		    }
		}
	    } else if (Mlsame(side, "R")) {
//Form C * H or C * H'  where  C = ( C1  C2 )
//W:= C * V'  =  (C1*V1' + C2 * V2')  (stored in WORK)
//W:= C2
		for (j = 0; j < k; j++) {
		    Ccopy(m, &C[(n - k + j) * ldc], 1, &work[j * ldwork], 1);
		}
//W: = W * V2'
		Ctrmm("Right", "Lower", "Conjugate transpose", "Unit", m, k, (COMPLEX) One, &V[(n - k) * ldv], ldv, work, ldwork);
		if (n > k) {
//W:= W + C1 * V1'
		    Cgemm("No transpose", "Conjugate transpose", m, k, n - k, (COMPLEX) One, C, ldc, V, ldv, (COMPLEX) One, work, ldwork);
		}
//W:= W * T or W * T'
		Ctrmm("Right", "Lower", trans, "Non-unit", m, k, (COMPLEX) One, T, ldt, work, ldwork);
//C:= C - W * V
		if (n > k) {
//C1:= C1 - W * V1
		    Cgemm("No transpose", "No transpose", m, n - k, k, (COMPLEX) - One, work, ldwork, V, ldv, (COMPLEX) One, C, ldc);
		}
//W:=W * V2
		Ctrmm("Right", "Lower", "No transpose", "Unit", m, k, (COMPLEX) One, &V[(n - k) * ldv], ldv, work, ldwork);
//C1: = C1 - W
		for (j = 0; j < k; j++) {
		    for (i = 0; i < m; i++) {
			C[i + (n - k + j) * ldc] = C[i + (n - k + j) * ldc] - work[i + j * ldwork];
		    }
		}
	    }
	}
    }
    return;
}
