/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlarfb.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rlarfb(const char *side, const char *trans, const char *direct,
       const char *storev, INTEGER m, INTEGER n, INTEGER k, REAL * V, INTEGER ldv, REAL * T, INTEGER ldt, REAL * C, INTEGER ldc, REAL * work, INTEGER ldwork)
{
    INTEGER i, j;
    REAL One = 1.0;
    char transt;

    //Quick return if possible
    if (m <= 0 || n <= 0) {
	return;
    }
    if (Mlsame(trans, "N")) {
	transt = 'T';
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
		    Rcopy(n, &C[j], ldc, &work[j * ldwork], 1);
		}
//W: = W * V1
		Rtrmm("Right", "Lower", "No transpose", "Unit", n, k, One, V, ldv, work, ldwork);
		if (m > k) {
//W: = W + C2'*V2
		    Rgemm("Transpose", "No transpose", n, k, m - k, One, &C[k], ldc, &V[k], ldv, One, work, ldwork);
		}
//W: = W * T'  or  W * T
		Rtrmm("Right", "Upper", &transt, "Non-unit", n, k, One, T, ldt, work, ldwork);
//C: = C - V * W'
		if (m > k) {
//C2: = C2 - V2 * W'
		    Rgemm("No transpose", "Transpose", m - k, n, k, -One, &V[k], ldv, work, ldwork, One, &C[k], ldc);
		}
//W: = W * V1'
		Rtrmm("Right", "Lower", "Transpose", "Unit", n, k, One, V, ldv, work, ldwork);
//C1: = C1 - W'
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[j + i * ldc] = C[j + i * ldc] - work[i + j * ldwork];
		    }
		}
	    } else if (Mlsame(side, "R")) {
//Form C * H or C * H'  where  C = ( C1  C2 )
//W: = C * V = (C1 * V1 + C2 * V2) (stored in WORK)
//W: = C1
		for (j = 0; j < k; j++) {
		    Rcopy(m, &C[j * ldc], 1, &work[j * ldwork], 1);
		}
//W: = W * V1
		Rtrmm("Right", "Lower", "No transpose", "Unit", m, k, One, V, ldv, work, ldwork);
		if (n > k) {
//W: = W + C2 * V2
		    Rgemm("No transpose", "No transpose", m, k, n - k, One, &C[k * ldc], ldc, &V[k], ldv, One, work, ldwork);
		}
//W: = W * T or W * T'
		Rtrmm("Right", "Upper", trans, "Non-unit", m, k, One, T, ldt, work, ldwork);
//C: = C - W * V'
		if (n > k) {
//C2: = C2 - W * V2'
		    Rgemm("No transpose", "Transpose", m, n - k, k, -One, work, ldwork, &V[k], ldv, One, &C[k * ldc], ldc);
		}
//W: = W * V1'
		Rtrmm("Right", "Lower", "Transpose", "Unit", m, k, One, V, ldv, work, ldwork);
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
		    Rcopy(n, &C[m - k + j], ldc, &work[j * ldwork], 1);
		}
//W: = W * V2
		Rtrmm("Right", "Upper", "No transpose", "Unit", n, k, One, &V[m - k], ldv, work, ldwork);
		if (m > k) {
//W: = W + C1'*V1
		    Rgemm("Transpose", "No transpose", n, k, m - k, One, C, ldc, V, ldv, One, work, ldwork);
		}
//W: = W * T'  or  W * T
		Rtrmm("Right", "Lower", &transt, "Non-unit", n, k, One, T, ldt, work, ldwork);
//C: = C - V * W'
		if (m > k) {
//C1:= C1 - V1 * W'
		    Rgemm("No transpose", "Transpose", m - k, n, k, -One, V, ldv, work, ldwork, One, C, ldc);
		}
//W: = W * V2'
		Rtrmm("Right", "Upper", "Transpose", "Unit", n, k, One, &V[m - k], ldv, work, ldwork);
//C2:= C2 - W'
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[m - k + j + i * ldc] = C[m - k + j + i * ldc] - work[i + j * ldwork];
		    }
		}
	    } else if (Mlsame(side, "R")) {
//Form C * H or C * H'  where  C = ( C1  C2 )
//W: = C * V = (C1 * V1 + C2 * V2) (stored in WORK)
//W: = C2
		for (j = 0; j < k; j++) {
		    Rcopy(m, &C[(n - k + j) * ldc], 1, &work[j * ldwork], 1);
		}
//W:= W * V2
		Rtrmm("Right", "Upper", "No transpose", "Unit", m, k, One, &V[n - k], ldv, work, ldwork);
		if (n > k) {
//W:= W + C1 * V1
		    Rgemm("No transpose", "No transpose", m, k, n - k, One, C, ldc, V, ldv, One, work, ldwork);
		}
//W:= W * T or W * T
		Rtrmm("Right", "Lower", trans, "Non-unit", m, k, One, T, ldt, work, ldwork);
//C:= C - W * V'
		if (n > k) {
//C1:= C1 - W * V1'
		    Rgemm("No transpose", "Transpose", m, n - k, k, -One, work, ldwork, V, ldv, One, C, ldc);
		}
//W: = W * V2'
		Rtrmm("Right", "Upper", "Transpose", "Unit", m, k, One, &V[n - k], ldv, work, ldwork);
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
		    Rcopy(n, &C[j], ldc, &work[j * ldwork], 1);
		}
//W:= W * V1'
		Rtrmm("Right", "Upper", "Transpose", "Unit", n, k, One, V, ldv, work, ldwork);
		if (m > k) {
//W:= W + C2'*V2'
		    Rgemm("Transpose", "Transpose", n, k, m - k, One, &C[k], ldc, &V[k * ldv], ldv, One, work, ldwork);
		}
//W:= W * T'  or  W * T
		Rtrmm("Right", "Upper", &transt, "Non-unit", n, k, One, T, ldt, work, ldwork);
//C:= C - V' * W'
		if (m > k) {
//C2:= C2 - V2' * W'
		    Rgemm("Transpose", "Transpose", m - k, n, k, -One, &V[k * ldv], ldv, work, ldwork, One, &C[k], ldc);
		}
//W:= W * V1
		Rtrmm("Right", "Upper", "No transpose", "Unit", n, k, One, V, ldv, work, ldwork);
//C1:= C1 - W'
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[j + i * ldc] = C[j + i * ldc] - work[i + j * ldwork];
		    }
		}
	    } else if (Mlsame(side, "R")) {
//Form C * H or C * H'  where  C = ( C1  C2 )
//W:= C * V'  =  (C1*V1' + C2 * V2')  (stored in WORK)
//W:= C1
		for (j = 0; j < k; j++) {
		    Rcopy(m, &C[j * ldc], 1, &work[j * ldwork], 1);
		}
//W:= W * V1'
		Rtrmm("Right", "Upper", "Transpose", "Unit", m, k, One, V, ldv, work, ldwork);
		if (n > k) {
//W:= W + C2 * V2'
		    Rgemm("No transpose", "Transpose", m, k, n - k, One, &C[k * ldc], ldc, &V[k * ldv], ldv, One, work, ldwork);
		}
//W:= W * T or W * T'
		Rtrmm("Right", "Upper", trans, "Non-unit", m, k, One, T, ldt, work, ldwork);
//C:= C - W * V
		if (n > k) {
//C2:= C2 - W * V2
		    Rgemm("No transpose", "No transpose", m, n - k, k, -One, work, ldwork, &V[k * ldv], ldv, One, &C[k * ldc], ldc);
		}
//W:= W * V1
		Rtrmm("Right", "Upper", "No transpose", "Unit", m, k, One, V, ldv, work, ldwork);
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
		    Rcopy(n, &C[m - k + j], ldc, &work[j * ldwork], 1);
		}
//W:= W * V2'
		Rtrmm("Right", "Lower", "Transpose", "Unit", n, k, One, &V[(m - k) * ldv], ldv, work, ldwork);
		if (m > k) {
//W:= W + C1'*V1'
		    Rgemm("Transpose", "Transpose", n, k, m - k, One, C, ldc, V, ldv, One, work, ldwork);
		}
//W:= W * T'  or  W * T
		Rtrmm("Right", "Lower", &transt, "Non-unit", n, k, One, T, ldt, work, ldwork);
//C:= C - V' * W'
		if (m > k) {
//C1:= C1 - V1' * W'
		    Rgemm("Transpose", "Transpose", m - k, n, k, -One, V, ldv, work, ldwork, One, C, ldc);
		}
//W:= W * V2
		Rtrmm("Right", "Lower", "No transpose", "Unit", n, k, One, &V[(m - k) * ldv], ldv, work, ldwork);
//C2:= C2 - W'
		for (j = 0; j < k; j++) {
		    for (i = 0; i < n; i++) {
			C[m - k + j + i * ldc] = C[m - k + j + i * ldc] - work[i + j * ldwork];
		    }
		}
	    } else if (Mlsame(side, "R")) {
//Form C * H or C * H'  where  C = ( C1  C2 )
//W:= C * V'  =  (C1*V1' + C2 * V2')  (stored in WORK)
//W:= C2
		for (j = 0; j < k; j++) {
		    Rcopy(m, &C[(n - k + j) * ldc], 1, &work[j * ldwork], 1);
		}
//W: = W * V2'
		Rtrmm("Right", "Lower", "Transpose", "Unit", m, k, One, &V[(n - k) * ldv], ldv, work, ldwork);
		if (n > k) {
//W:= W + C1 * V1'
		    Rgemm("No transpose", "Transpose", m, k, n - k, One, C, ldc, V, ldv, One, work, ldwork);
		}
//W:= W * T or W * T'
		Rtrmm("Right", "Lower", trans, "Non-unit", m, k, One, T, ldt, work, ldwork);
//C:= C - W * V
		if (n > k) {
//C1:= C1 - W * V1
		    Rgemm("No transpose", "No transpose", m, n - k, k, -One, work, ldwork, V, ldv, One, C, ldc);
		}
//W:=W * V2
		Rtrmm("Right", "Lower", "No transpose", "Unit", m, k, One, &V[(n - k) * ldv], ldv, work, ldwork);
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
