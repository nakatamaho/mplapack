/*
 * Copyright (c) 2010-2011
 *      RIKEN
 *      All rights reserved.
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
  Contributed by Takao, Yasuyoshi and Nakata, Maho, 2010-2011
*/


#ifdef __DD_REAL_CUDA_H__
#else
#define __DD_REAL_CUDA_H__
#include  <cstdlib>
struct dd_real {
  double x[2];
};

#define _QD_SPLITTER 134217729.0 // 2^27 + 1 
#define _SD_SPLITTER 4096.0 // 2^12 + 1 

// Branch-free dd arithmetic of libQD3 built with QD_BF: dd_add is
// dd_real::bf_add (Kouya Algorithm 6, DWBFAdd), dd_mul dd_real::bf_mul
// (Algorithm 8, DWBFMul), with the error-free product computed by fma.
// Compile with -fmad=false.

inline __host__ __device__ double dd_two_sum(double a, double b, double &err)
{
    double s = a + b;
    double bb = s - a;
    err = (a - (s - bb)) + (b - bb);
    return s;
}

inline __host__ __device__ double dd_quick_two_sum(double a, double b, double &err)
{
    double s = a + b;
    err = b - (s - a);
    return s;
}

inline __host__ __device__ void dd_bf_mul(const dd_real &a, const dd_real &b, dd_real &c)
{
    double p00, pe00, p01, p10, g1, g2;
    p00 = a.x[0] * b.x[0];
    pe00 = fma(a.x[0], b.x[0], -p00);
    p01 = a.x[0] * b.x[1];
    p10 = a.x[1] * b.x[0];
    g1 = p01 + p10;
    g2 = pe00 + g1;
    c.x[0] = dd_quick_two_sum(p00, g2, c.x[1]);
}

inline __host__ __device__ void dd_bf_add(const dd_real &a, const dd_real &b, dd_real &c)
{
    double g1, g1e, g2, g2e, g3, g3e, g4, g5;
    g1 = dd_two_sum(a.x[0], b.x[0], g1e);
    g2 = dd_two_sum(a.x[1], b.x[1], g2e);
    g3 = dd_quick_two_sum(g1, g2, g3e);
    g4 = g1e + g2e;
    g5 = g4 + g3e;
    c.x[0] = dd_quick_two_sum(g3, g5, c.x[1]);
}

inline __device__ void dd_mul(const dd_real &a, const dd_real &b, dd_real &c) { dd_bf_mul(a, b, c); }
inline __device__ void dd_add(const dd_real &a, const dd_real &b, dd_real &c) { dd_bf_add(a, b, c); }
inline void dd_mul_host(const dd_real &a, const dd_real &b, dd_real &c) { dd_bf_mul(a, b, c); }
inline void dd_add_host(const dd_real &a, const dd_real &b, dd_real &c) { dd_bf_add(a, b, c); }

inline __device__ void dd_mad(dd_real &c, const dd_real &a, const dd_real &b)
{
    dd_real temp;
    
    dd_mul( a, b, temp);
    dd_add( temp, c, c);
}

inline void dd_set(dd_real &c, const double a, const double b)
{
    c.x[0] = a;
    c.x[1] = b;
}

inline int dd_eq(const dd_real &c1, const dd_real &c2)
{
    return c1.x[0] == c2.x[0] && c1.x[1] == c2.x[1];
}

#endif
