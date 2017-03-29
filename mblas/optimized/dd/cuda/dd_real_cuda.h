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

#define SLOPPY_MUL
#define CRAY_ADD

#define _QD_SPLITTER 134217729.0 // 2^27 + 1 
#define _SD_SPLITTER 4096.0 // 2^12 + 1 

// SLOPPY_MUL and CRAY_ADD are available
// therefore, all possible combinations are four.

#ifndef SLOPPY_MUL
inline __device__ void dd_mul(const dd_real &a, const dd_real &b, dd_real &c)
{
    double p;
    double t;
    double a_hh,b_hh;
    double a_hl,b_hl;
    double e;

    p = __dmul_rn(a.x[0], b.x[0]);
    e = __fma_rn(a.x[0], b.x[0], -1.0 * p);
    e = e + __fma_rn ( a.x[1], b.x[0],  __dmul_rn(a.x[0],b.x[1]));

    c.x[0] = p + e;
    c.x[1] = e - (c.x[0] - p);
}

#else

inline __device__ void dd_mul(const dd_real &a, const dd_real &b, dd_real &c)
{
    double m1,m2,t,ph,e,pl;

    m1 = a.x[0] * b.x[1];
    m2 = a.x[1] * b.x[0];
    
    t = m1 + m2;
    ph = __fma_rn(a.x[0], b.x[0], t);
    e = __fma_rn(a.x[0], b.x[0], -1.0 * ph); 
    pl = e + t;
    
    c.x[0] = ph;
    c.x[1] = pl;
}

#endif

inline void dd_mul_host(const dd_real &a, const dd_real &b, dd_real &c)
{
    double p;
    double t;
    double a_hh,b_hh;
    double a_hl,b_hl;
    double e;

    p = a.x[0] * b.x[0];
    t = 134217729.0 * a.x[0];
    a_hh = t - (t - a.x[0]);
    a_hl = a.x[0] - a_hh;
    t = 134217729.0 * b.x[0];
    b_hh = t - (t - b.x[0]);
    b_hl = b.x[0] - b_hh;
    e = ((a_hh * b_hh - p) + a_hh * b_hl + a_hl * b_hh) + a_hl * b_hl;
    e = e + (a.x[0] * b.x[1] + a.x[1] * b.x[0]);
    
    c.x[0] = p + e;
    c.x[1] = e - (c.x[0] - p);
}

#ifndef CRAY_ADD

inline __device__ void dd_add(const dd_real &a, const dd_real &b, dd_real &c)
{
    double s1,s2;
    double v;
    double t1,t2;

    s1 = a.x[0] + b.x[0];
    v = s1 - a.x[0];
    s2 = ((b.x[0] - v) + (a.x[0] - (s1 - v)));
    t1 = a.x[1] + b.x[1];
    v = t1 - a.x[1];
    t2 = ((b.x[1] - v) + (a.x[1] - (t1 - v)));
    s2 = s2 + t1;
    t1 = s1 + s2;
    s2 = s2 - (t1 - s1);
    t2 = t2 + s2;
    
    c.x[0] = t1 + t2;
    c.x[1] = t2 - (c.x[0] - t1);
        
}

#else

inline __device__ void dd_add(const dd_real &a, const dd_real &b, dd_real &c)
{
    double s1,s2;
    double v;
    double t1,t2;

    //s = qd::two_sum(a.x[0], b.x[0], e);
    s1 = a.x[0] + b.x[0];
    v = s1 - a.x[0];
    s2 = ((b.x[0] - v) + (a.x[0] - (s1 - v)));

    //e += (a.x[1] + b.x[1]);
    s2 = s2 + (a.x[1] + b.x[1]);

    //s = qd::quick_two_sum(s, e, e);
    //return dd_real(s, e);
    t1 = s1 + s2;
    t2 = s2 - (t1 - s1);

    c.x[0] = t1;
    c.x[1] = t2;
}

#endif

inline void dd_add_host(const dd_real &a, const dd_real &b, dd_real &c)
{
    double s1,s2;
    double v;
    double t1,t2;

    s1 = a.x[0] + b.x[0];
    v = s1 - a.x[0];
    s2 = ((b.x[0] - v) + (a.x[0] - (s1 - v)));
    t1 = a.x[1] + b.x[1];
    v = t1 - a.x[1];
    t2 = ((b.x[1] - v) + (a.x[1] - (t1 - v)));
    s2 = s2 + t1;
    t1 = s1 + s2;
    s2 = s2 - (t1 - s1);
    t2 = t2 + s2;
    
    c.x[0] = t1 + t2;
    c.x[1] = t2 - (c.x[0] - t1);
}

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
