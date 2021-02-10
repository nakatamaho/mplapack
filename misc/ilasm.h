/*
 * Copyright (c) 2010
 *	Yuichiro Yasui <yasui@indsys.chuo-u.ac.jp>
 * 	All rights reserved.
 *
 * $Id: ilasm.h,v 1.1 2011/01/05 07:32:51 nakatamaho Exp $
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

#ifndef __inline_asm_h__
#define __inline_asm_h__


#if defined (__cplusplus)
extern "C" {
#endif

  /* ------------------------------------------------------------
   * rdtsc: Read Time-Stamp Counter
   * 
   * USAGE:
   *   unsigned long t1, t2, ts;
   *   t1 = rdtsc();
   *   target();
   *   t2 = rdtsc();
   *   ts = t2 - t1;
   * ------------------------------------------------------------ */
  static __inline__ unsigned long rdtsc(void) {
    unsigned long a, d;
    __asm__ __volatile__ ("rdtsc" : "=a" (a), "=d" (d));
    return a | (d << 32);
  };
  
  /* ------------------------------------------------------------
   * apicid: /proc/meminfo
   * 
   * USAGE:
   *   int id;
   *   id = get_apicid();
   * ------------------------------------------------------------ */
  static __inline__ int get_apicid(void) {
    int eax, ebx, ecx, edx;
    __asm__ __volatile__ ("cpuid" : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) : "0" (1));
    return ((ebx >> 24) & 0xff);
  }

#if defined (__cplusplus)
}
#endif

#endif	/* __inline_asm_h__ */
