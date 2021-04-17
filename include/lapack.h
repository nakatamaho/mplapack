/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: lapack.h,v 1.9 2010/08/07 03:15:46 nakatamaho Exp $
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
  LAPACK header prototypes based on version 3.1.1 are here
  Define to a macro mangling the given C identifier (in lower and upper
  case), which must not contain underscores, for linking with Fortran.
*/

#ifndef _LAPACK_H_
#define _LAPACK_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
typedef logical (*L_fp)(...);
#else
typedef logical(*L_fp);
#endif

/* LAPACK headers (3.1.1) */
#define cbdsqr_f77 F77_FUNC(cbdsqr, CBDSQR)
#define cgbbrd_f77 F77_FUNC(cgbbrd, CGBBRD)
#define cgbcon_f77 F77_FUNC(cgbcon, CGBCON)
#define cgbequ_f77 F77_FUNC(cgbequ, CGBEQU)
#define cgbrfs_f77 F77_FUNC(cgbrfs, CGBRFS)
#define cgbsv_f77 F77_FUNC(cgbsv, CGBSV)
#define cgbsvx_f77 F77_FUNC(cgbsvx, CGBSVX)
#define cgbtf2_f77 F77_FUNC(cgbtf2, CGBTF2)
#define cgbtrf_f77 F77_FUNC(cgbtrf, CGBTRF)
#define cgbtrs_f77 F77_FUNC(cgbtrs, CGBTRS)
#define cgebak_f77 F77_FUNC(cgebak, CGEBAK)
#define cgebal_f77 F77_FUNC(cgebal, CGEBAL)
#define cgebd2_f77 F77_FUNC(cgebd2, CGEBD2)
#define cgebrd_f77 F77_FUNC(cgebrd, CGEBRD)
#define cgecon_f77 F77_FUNC(cgecon, CGECON)
#define cgeequ_f77 F77_FUNC(cgeequ, CGEEQU)
#define cgees_f77 F77_FUNC(cgees, CGEES)
#define cgeesx_f77 F77_FUNC(cgeesx, CGEESX)
#define cgeev_f77 F77_FUNC(cgeev, CGEEV)
#define cgeevx_f77 F77_FUNC(cgeevx, CGEEVX)
#define cgegs_f77 F77_FUNC(cgegs, CGEGS)
#define cgegv_f77 F77_FUNC(cgegv, CGEGV)
#define cgehd2_f77 F77_FUNC(cgehd2, CGEHD2)
#define cgehrd_f77 F77_FUNC(cgehrd, CGEHRD)
#define cgelq2_f77 F77_FUNC(cgelq2, CGELQ2)
#define cgelqf_f77 F77_FUNC(cgelqf, CGELQF)
#define cgels_f77 F77_FUNC(cgels, CGELS)
#define cgelsd_f77 F77_FUNC(cgelsd, CGELSD)
#define cgelss_f77 F77_FUNC(cgelss, CGELSS)
#define cgelsx_f77 F77_FUNC(cgelsx, CGELSX)
#define cgelsy_f77 F77_FUNC(cgelsy, CGELSY)
#define cgeql2_f77 F77_FUNC(cgeql2, CGEQL2)
#define cgeqlf_f77 F77_FUNC(cgeqlf, CGEQLF)
#define cgeqp3_f77 F77_FUNC(cgeqp3, CGEQP3)
#define cgeqpf_f77 F77_FUNC(cgeqpf, CGEQPF)
#define cgeqr2_f77 F77_FUNC(cgeqr2, CGEQR2)
#define cgeqrf_f77 F77_FUNC(cgeqrf, CGEQRF)
#define cgerfs_f77 F77_FUNC(cgerfs, CGERFS)
#define cgerq2_f77 F77_FUNC(cgerq2, CGERQ2)
#define cgerqf_f77 F77_FUNC(cgerqf, CGERQF)
#define cgesc2_f77 F77_FUNC(cgesc2, CGESC2)
#define cgesdd_f77 F77_FUNC(cgesdd, CGESDD)
#define cgesv_f77 F77_FUNC(cgesv, CGESV)
#define cgesvd_f77 F77_FUNC(cgesvd, CGESVD)
#define cgesvx_f77 F77_FUNC(cgesvx, CGESVX)
#define cgetc2_f77 F77_FUNC(cgetc2, CGETC2)
#define cgetf2_f77 F77_FUNC(cgetf2, CGETF2)
#define cgetrf_f77 F77_FUNC(cgetrf, CGETRF)
#define cgetri_f77 F77_FUNC(cgetri, CGETRI)
#define cgetrs_f77 F77_FUNC(cgetrs, CGETRS)
#define cggbak_f77 F77_FUNC(cggbak, CGGBAK)
#define cggbal_f77 F77_FUNC(cggbal, CGGBAL)
#define cgges_f77 F77_FUNC(cgges, CGGES)
#define cggesx_f77 F77_FUNC(cggesx, CGGESX)
#define cggev_f77 F77_FUNC(cggev, CGGEV)
#define cggevx_f77 F77_FUNC(cggevx, CGGEVX)
#define cggglm_f77 F77_FUNC(cggglm, CGGGLM)
#define cgghrd_f77 F77_FUNC(cgghrd, CGGHRD)
#define cgglse_f77 F77_FUNC(cgglse, CGGLSE)
#define cggqrf_f77 F77_FUNC(cggqrf, CGGQRF)
#define cggrqf_f77 F77_FUNC(cggrqf, CGGRQF)
#define cggsvd_f77 F77_FUNC(cggsvd, CGGSVD)
#define cggsvp_f77 F77_FUNC(cggsvp, CGGSVP)
#define cgtcon_f77 F77_FUNC(cgtcon, CGTCON)
#define cgtrfs_f77 F77_FUNC(cgtrfs, CGTRFS)
#define cgtsv_f77 F77_FUNC(cgtsv, CGTSV)
#define cgtsvx_f77 F77_FUNC(cgtsvx, CGTSVX)
#define cgttrf_f77 F77_FUNC(cgttrf, CGTTRF)
#define cgttrs_f77 F77_FUNC(cgttrs, CGTTRS)
#define cgtts2_f77 F77_FUNC(cgtts2, CGTTS2)
#define chbev_f77 F77_FUNC(chbev, CHBEV)
#define chbevd_f77 F77_FUNC(chbevd, CHBEVD)
#define chbevx_f77 F77_FUNC(chbevx, CHBEVX)
#define chbgst_f77 F77_FUNC(chbgst, CHBGST)
#define chbgv_f77 F77_FUNC(chbgv, CHBGV)
#define chbgvd_f77 F77_FUNC(chbgvd, CHBGVD)
#define chbgvx_f77 F77_FUNC(chbgvx, CHBGVX)
#define chbtrd_f77 F77_FUNC(chbtrd, CHBTRD)
#define checon_f77 F77_FUNC(checon, CHECON)
#define cheev_f77 F77_FUNC(cheev, CHEEV)
#define cheevd_f77 F77_FUNC(cheevd, CHEEVD)
#define cheevr_f77 F77_FUNC(cheevr, CHEEVR)
#define cheevx_f77 F77_FUNC(cheevx, CHEEVX)
#define chegs2_f77 F77_FUNC(chegs2, CHEGS2)
#define chegst_f77 F77_FUNC(chegst, CHEGST)
#define chegv_f77 F77_FUNC(chegv, CHEGV)
#define chegvd_f77 F77_FUNC(chegvd, CHEGVD)
#define chegvx_f77 F77_FUNC(chegvx, CHEGVX)
#define cherfs_f77 F77_FUNC(cherfs, CHERFS)
#define chesv_f77 F77_FUNC(chesv, CHESV)
#define chesvx_f77 F77_FUNC(chesvx, CHESVX)
#define chetd2_f77 F77_FUNC(chetd2, CHETD2)
#define chetf2_f77 F77_FUNC(chetf2, CHETF2)
#define chetrd_f77 F77_FUNC(chetrd, CHETRD)
#define chetrf_f77 F77_FUNC(chetrf, CHETRF)
#define chetri_f77 F77_FUNC(chetri, CHETRI)
#define chetrs_f77 F77_FUNC(chetrs, CHETRS)
#define chgeqz_f77 F77_FUNC(chgeqz, CHGEQZ)
#define chpcon_f77 F77_FUNC(chpcon, CHPCON)
#define chpev_f77 F77_FUNC(chpev, CHPEV)
#define chpevd_f77 F77_FUNC(chpevd, CHPEVD)
#define chpevx_f77 F77_FUNC(chpevx, CHPEVX)
#define chpgst_f77 F77_FUNC(chpgst, CHPGST)
#define chpgv_f77 F77_FUNC(chpgv, CHPGV)
#define chpgvd_f77 F77_FUNC(chpgvd, CHPGVD)
#define chpgvx_f77 F77_FUNC(chpgvx, CHPGVX)
#define chprfs_f77 F77_FUNC(chprfs, CHPRFS)
#define chpsv_f77 F77_FUNC(chpsv, CHPSV)
#define chpsvx_f77 F77_FUNC(chpsvx, CHPSVX)
#define chptrd_f77 F77_FUNC(chptrd, CHPTRD)
#define chptrf_f77 F77_FUNC(chptrf, CHPTRF)
#define chptri_f77 F77_FUNC(chptri, CHPTRI)
#define chptrs_f77 F77_FUNC(chptrs, CHPTRS)
#define chsein_f77 F77_FUNC(chsein, CHSEIN)
#define chseqr_f77 F77_FUNC(chseqr, CHSEQR)
#define clabrd_f77 F77_FUNC(clabrd, CLABRD)
#define clacgv_f77 F77_FUNC(clacgv, CLACGV)
#define clacn2_f77 F77_FUNC(clacn2, CLACN2)
#define clacon_f77 F77_FUNC(clacon, CLACON)
#define clacp2_f77 F77_FUNC(clacp2, CLACP2)
#define clacpy_f77 F77_FUNC(clacpy, CLACPY)
#define clacrm_f77 F77_FUNC(clacrm, CLACRM)
#define clacrt_f77 F77_FUNC(clacrt, CLACRT)
#define cladiv_f77 F77_FUNC(cladiv, CLADIV)
#define claed0_f77 F77_FUNC(claed0, CLAED0)
#define claed7_f77 F77_FUNC(claed7, CLAED7)
#define claed8_f77 F77_FUNC(claed8, CLAED8)
#define claein_f77 F77_FUNC(claein, CLAEIN)
#define claesy_f77 F77_FUNC(claesy, CLAESY)
#define claev2_f77 F77_FUNC(claev2, CLAEV2)
#define clag2z_f77 F77_FUNC(clag2z, CLAG2Z)
#define clags2_f77 F77_FUNC(clags2, CLAGS2)
#define clagtm_f77 F77_FUNC(clagtm, CLAGTM)
#define clahef_f77 F77_FUNC(clahef, CLAHEF)
#define clahqr_f77 F77_FUNC(clahqr, CLAHQR)
#define clahr2_f77 F77_FUNC(clahr2, CLAHR2)
#define clahrd_f77 F77_FUNC(clahrd, CLAHRD)
#define claic1_f77 F77_FUNC(claic1, CLAIC1)
#define clals0_f77 F77_FUNC(clals0, CLALS0)
#define clalsa_f77 F77_FUNC(clalsa, CLALSA)
#define clalsd_f77 F77_FUNC(clalsd, CLALSD)
#define clangb_f77 F77_FUNC(clangb, CLANGB)
#define clange_f77 F77_FUNC(clange, CLANGE)
#define clangt_f77 F77_FUNC(clangt, CLANGT)
#define clanhb_f77 F77_FUNC(clanhb, CLANHB)
#define clanhe_f77 F77_FUNC(clanhe, CLANHE)
#define clanhp_f77 F77_FUNC(clanhp, CLANHP)
#define clanhs_f77 F77_FUNC(clanhs, CLANHS)
#define clanht_f77 F77_FUNC(clanht, CLANHT)
#define clansb_f77 F77_FUNC(clansb, CLANSB)
#define clansp_f77 F77_FUNC(clansp, CLANSP)
#define clansy_f77 F77_FUNC(clansy, CLANSY)
#define clantb_f77 F77_FUNC(clantb, CLANTB)
#define clantp_f77 F77_FUNC(clantp, CLANTP)
#define clantr_f77 F77_FUNC(clantr, CLANTR)
#define clapll_f77 F77_FUNC(clapll, CLAPLL)
#define clapmt_f77 F77_FUNC(clapmt, CLAPMT)
#define claqgb_f77 F77_FUNC(claqgb, CLAQGB)
#define claqge_f77 F77_FUNC(claqge, CLAQGE)
#define claqhb_f77 F77_FUNC(claqhb, CLAQHB)
#define claqhe_f77 F77_FUNC(claqhe, CLAQHE)
#define claqhp_f77 F77_FUNC(claqhp, CLAQHP)
#define claqp2_f77 F77_FUNC(claqp2, CLAQP2)
#define claqps_f77 F77_FUNC(claqps, CLAQPS)
#define claqr0_f77 F77_FUNC(claqr0, CLAQR0)
#define claqr1_f77 F77_FUNC(claqr1, CLAQR1)
#define claqr2_f77 F77_FUNC(claqr2, CLAQR2)
#define claqr3_f77 F77_FUNC(claqr3, CLAQR3)
#define claqr4_f77 F77_FUNC(claqr4, CLAQR4)
#define claqr5_f77 F77_FUNC(claqr5, CLAQR5)
#define claqsb_f77 F77_FUNC(claqsb, CLAQSB)
#define claqsp_f77 F77_FUNC(claqsp, CLAQSP)
#define claqsy_f77 F77_FUNC(claqsy, CLAQSY)
#define clar1v_f77 F77_FUNC(clar1v, CLAR1V)
#define clar2v_f77 F77_FUNC(clar2v, CLAR2V)
#define clarcm_f77 F77_FUNC(clarcm, CLARCM)
#define clarf_f77 F77_FUNC(clarf, CLARF)
#define clarfb_f77 F77_FUNC(clarfb, CLARFB)
#define clarfg_f77 F77_FUNC(clarfg, CLARFG)
#define clarft_f77 F77_FUNC(clarft, CLARFT)
#define clarfx_f77 F77_FUNC(clarfx, CLARFX)
#define clargv_f77 F77_FUNC(clargv, CLARGV)
#define clarnv_f77 F77_FUNC(clarnv, CLARNV)
#define clarrv_f77 F77_FUNC(clarrv, CLARRV)
#define clartg_f77 F77_FUNC(clartg, CLARTG)
#define clartv_f77 F77_FUNC(clartv, CLARTV)
#define clarz_f77 F77_FUNC(clarz, CLARZ)
#define clarzb_f77 F77_FUNC(clarzb, CLARZB)
#define clarzt_f77 F77_FUNC(clarzt, CLARZT)
#define clascl_f77 F77_FUNC(clascl, CLASCL)
#define claset_f77 F77_FUNC(claset, CLASET)
#define clasr_f77 F77_FUNC(clasr, CLASR)
#define classq_f77 F77_FUNC(classq, CLASSQ)
#define claswp_f77 F77_FUNC(claswp, CLASWP)
#define clasyf_f77 F77_FUNC(clasyf, CLASYF)
#define clatbs_f77 F77_FUNC(clatbs, CLATBS)
#define clatdf_f77 F77_FUNC(clatdf, CLATDF)
#define clatps_f77 F77_FUNC(clatps, CLATPS)
#define clatrd_f77 F77_FUNC(clatrd, CLATRD)
#define clatrs_f77 F77_FUNC(clatrs, CLATRS)
#define clatrz_f77 F77_FUNC(clatrz, CLATRZ)
#define clatzm_f77 F77_FUNC(clatzm, CLATZM)
#define clauu2_f77 F77_FUNC(clauu2, CLAUU2)
#define clauum_f77 F77_FUNC(clauum, CLAUUM)
#define cpbcon_f77 F77_FUNC(cpbcon, CPBCON)
#define cpbequ_f77 F77_FUNC(cpbequ, CPBEQU)
#define cpbrfs_f77 F77_FUNC(cpbrfs, CPBRFS)
#define cpbstf_f77 F77_FUNC(cpbstf, CPBSTF)
#define cpbsv_f77 F77_FUNC(cpbsv, CPBSV)
#define cpbsvx_f77 F77_FUNC(cpbsvx, CPBSVX)
#define cpbtf2_f77 F77_FUNC(cpbtf2, CPBTF2)
#define cpbtrf_f77 F77_FUNC(cpbtrf, CPBTRF)
#define cpbtrs_f77 F77_FUNC(cpbtrs, CPBTRS)
#define cpocon_f77 F77_FUNC(cpocon, CPOCON)
#define cpoequ_f77 F77_FUNC(cpoequ, CPOEQU)
#define cporfs_f77 F77_FUNC(cporfs, CPORFS)
#define cposv_f77 F77_FUNC(cposv, CPOSV)
#define cposvx_f77 F77_FUNC(cposvx, CPOSVX)
#define cpotf2_f77 F77_FUNC(cpotf2, CPOTF2)
#define cpotrf_f77 F77_FUNC(cpotrf, CPOTRF)
#define cpotri_f77 F77_FUNC(cpotri, CPOTRI)
#define cpotrs_f77 F77_FUNC(cpotrs, CPOTRS)
#define cppcon_f77 F77_FUNC(cppcon, CPPCON)
#define cppequ_f77 F77_FUNC(cppequ, CPPEQU)
#define cpprfs_f77 F77_FUNC(cpprfs, CPPRFS)
#define cppsv_f77 F77_FUNC(cppsv, CPPSV)
#define cppsvx_f77 F77_FUNC(cppsvx, CPPSVX)
#define cpptrf_f77 F77_FUNC(cpptrf, CPPTRF)
#define cpptri_f77 F77_FUNC(cpptri, CPPTRI)
#define cpptrs_f77 F77_FUNC(cpptrs, CPPTRS)
#define cptcon_f77 F77_FUNC(cptcon, CPTCON)
#define cpteqr_f77 F77_FUNC(cpteqr, CPTEQR)
#define cptrfs_f77 F77_FUNC(cptrfs, CPTRFS)
#define cptsv_f77 F77_FUNC(cptsv, CPTSV)
#define cptsvx_f77 F77_FUNC(cptsvx, CPTSVX)
#define cpttrf_f77 F77_FUNC(cpttrf, CPTTRF)
#define cpttrs_f77 F77_FUNC(cpttrs, CPTTRS)
#define cptts2_f77 F77_FUNC(cptts2, CPTTS2)
#define crot_f77 F77_FUNC(crot, CROT)
#define cspcon_f77 F77_FUNC(cspcon, CSPCON)
#define cspmv_f77 F77_FUNC(cspmv, CSPMV)
#define cspr_f77 F77_FUNC(cspr, CSPR)
#define csprfs_f77 F77_FUNC(csprfs, CSPRFS)
#define cspsv_f77 F77_FUNC(cspsv, CSPSV)
#define cspsvx_f77 F77_FUNC(cspsvx, CSPSVX)
#define csptrf_f77 F77_FUNC(csptrf, CSPTRF)
#define csptri_f77 F77_FUNC(csptri, CSPTRI)
#define csptrs_f77 F77_FUNC(csptrs, CSPTRS)
#define csrscl_f77 F77_FUNC(csrscl, CSRSCL)
#define cstedc_f77 F77_FUNC(cstedc, CSTEDC)
#define cstegr_f77 F77_FUNC(cstegr, CSTEGR)
#define cstein_f77 F77_FUNC(cstein, CSTEIN)
#define cstemr_f77 F77_FUNC(cstemr, CSTEMR)
#define csteqr_f77 F77_FUNC(csteqr, CSTEQR)
#define csycon_f77 F77_FUNC(csycon, CSYCON)
#define csymv_f77 F77_FUNC(csymv, CSYMV)
#define csyr_f77 F77_FUNC(csyr, CSYR)
#define csyrfs_f77 F77_FUNC(csyrfs, CSYRFS)
#define csysv_f77 F77_FUNC(csysv, CSYSV)
#define csysvx_f77 F77_FUNC(csysvx, CSYSVX)
#define csytf2_f77 F77_FUNC(csytf2, CSYTF2)
#define csytrf_f77 F77_FUNC(csytrf, CSYTRF)
#define csytri_f77 F77_FUNC(csytri, CSYTRI)
#define csytrs_f77 F77_FUNC(csytrs, CSYTRS)
#define ctbcon_f77 F77_FUNC(ctbcon, CTBCON)
#define ctbrfs_f77 F77_FUNC(ctbrfs, CTBRFS)
#define ctbtrs_f77 F77_FUNC(ctbtrs, CTBTRS)
#define ctgevc_f77 F77_FUNC(ctgevc, CTGEVC)
#define ctgex2_f77 F77_FUNC(ctgex2, CTGEX2)
#define ctgexc_f77 F77_FUNC(ctgexc, CTGEXC)
#define ctgsen_f77 F77_FUNC(ctgsen, CTGSEN)
#define ctgsja_f77 F77_FUNC(ctgsja, CTGSJA)
#define ctgsna_f77 F77_FUNC(ctgsna, CTGSNA)
#define ctgsy2_f77 F77_FUNC(ctgsy2, CTGSY2)
#define ctgsyl_f77 F77_FUNC(ctgsyl, CTGSYL)
#define ctpcon_f77 F77_FUNC(ctpcon, CTPCON)
#define ctprfs_f77 F77_FUNC(ctprfs, CTPRFS)
#define ctptri_f77 F77_FUNC(ctptri, CTPTRI)
#define ctptrs_f77 F77_FUNC(ctptrs, CTPTRS)
#define ctrcon_f77 F77_FUNC(ctrcon, CTRCON)
#define ctrevc_f77 F77_FUNC(ctrevc, CTREVC)
#define ctrexc_f77 F77_FUNC(ctrexc, CTREXC)
#define ctrrfs_f77 F77_FUNC(ctrrfs, CTRRFS)
#define ctrsen_f77 F77_FUNC(ctrsen, CTRSEN)
#define ctrsna_f77 F77_FUNC(ctrsna, CTRSNA)
#define ctrsyl_f77 F77_FUNC(ctrsyl, CTRSYL)
#define ctrti2_f77 F77_FUNC(ctrti2, CTRTI2)
#define ctrtri_f77 F77_FUNC(ctrtri, CTRTRI)
#define ctrtrs_f77 F77_FUNC(ctrtrs, CTRTRS)
#define ctzrqf_f77 F77_FUNC(ctzrqf, CTZRQF)
#define ctzrzf_f77 F77_FUNC(ctzrzf, CTZRZF)
#define cung2l_f77 F77_FUNC(cung2l, CUNG2L)
#define cung2r_f77 F77_FUNC(cung2r, CUNG2R)
#define cungbr_f77 F77_FUNC(cungbr, CUNGBR)
#define cunghr_f77 F77_FUNC(cunghr, CUNGHR)
#define cungl2_f77 F77_FUNC(cungl2, CUNGL2)
#define cunglq_f77 F77_FUNC(cunglq, CUNGLQ)
#define cungql_f77 F77_FUNC(cungql, CUNGQL)
#define cungqr_f77 F77_FUNC(cungqr, CUNGQR)
#define cungr2_f77 F77_FUNC(cungr2, CUNGR2)
#define cungrq_f77 F77_FUNC(cungrq, CUNGRQ)
#define cungtr_f77 F77_FUNC(cungtr, CUNGTR)
#define cunm2l_f77 F77_FUNC(cunm2l, CUNM2L)
#define cunm2r_f77 F77_FUNC(cunm2r, CUNM2R)
#define cunmbr_f77 F77_FUNC(cunmbr, CUNMBR)
#define cunmhr_f77 F77_FUNC(cunmhr, CUNMHR)
#define cunml2_f77 F77_FUNC(cunml2, CUNML2)
#define cunmlq_f77 F77_FUNC(cunmlq, CUNMLQ)
#define cunmql_f77 F77_FUNC(cunmql, CUNMQL)
#define cunmqr_f77 F77_FUNC(cunmqr, CUNMQR)
#define cunmr2_f77 F77_FUNC(cunmr2, CUNMR2)
#define cunmr3_f77 F77_FUNC(cunmr3, CUNMR3)
#define cunmrq_f77 F77_FUNC(cunmrq, CUNMRQ)
#define cunmrz_f77 F77_FUNC(cunmrz, CUNMRZ)
#define cunmtr_f77 F77_FUNC(cunmtr, CUNMTR)
#define cupgtr_f77 F77_FUNC(cupgtr, CUPGTR)
#define cupmtr_f77 F77_FUNC(cupmtr, CUPMTR)
#define dbdsdc_f77 F77_FUNC(dbdsdc, DBDSDC)
#define dbdsqr_f77 F77_FUNC(dbdsqr, DBDSQR)
#define ddisna_f77 F77_FUNC(ddisna, DDISNA)
#define dgbbrd_f77 F77_FUNC(dgbbrd, DGBBRD)
#define dgbcon_f77 F77_FUNC(dgbcon, DGBCON)
#define dgbequ_f77 F77_FUNC(dgbequ, DGBEQU)
#define dgbrfs_f77 F77_FUNC(dgbrfs, DGBRFS)
#define dgbsv_f77 F77_FUNC(dgbsv, DGBSV)
#define dgbsvx_f77 F77_FUNC(dgbsvx, DGBSVX)
#define dgbtf2_f77 F77_FUNC(dgbtf2, DGBTF2)
#define dgbtrf_f77 F77_FUNC(dgbtrf, DGBTRF)
#define dgbtrs_f77 F77_FUNC(dgbtrs, DGBTRS)
#define dgebak_f77 F77_FUNC(dgebak, DGEBAK)
#define dgebal_f77 F77_FUNC(dgebal, DGEBAL)
#define dgebd2_f77 F77_FUNC(dgebd2, DGEBD2)
#define dgebrd_f77 F77_FUNC(dgebrd, DGEBRD)
#define dgecon_f77 F77_FUNC(dgecon, DGECON)
#define dgeequ_f77 F77_FUNC(dgeequ, DGEEQU)
#define dgees_f77 F77_FUNC(dgees, DGEES)
#define dgeesx_f77 F77_FUNC(dgeesx, DGEESX)
#define dgeev_f77 F77_FUNC(dgeev, DGEEV)
#define dgeevx_f77 F77_FUNC(dgeevx, DGEEVX)
#define dgegs_f77 F77_FUNC(dgegs, DGEGS)
#define dgegv_f77 F77_FUNC(dgegv, DGEGV)
#define dgehd2_f77 F77_FUNC(dgehd2, DGEHD2)
#define dgehrd_f77 F77_FUNC(dgehrd, DGEHRD)
#define dgelq2_f77 F77_FUNC(dgelq2, DGELQ2)
#define dgelqf_f77 F77_FUNC(dgelqf, DGELQF)
#define dgels_f77 F77_FUNC(dgels, DGELS)
#define dgelsd_f77 F77_FUNC(dgelsd, DGELSD)
#define dgelss_f77 F77_FUNC(dgelss, DGELSS)
#define dgelsx_f77 F77_FUNC(dgelsx, DGELSX)
#define dgelsy_f77 F77_FUNC(dgelsy, DGELSY)
#define dgeql2_f77 F77_FUNC(dgeql2, DGEQL2)
#define dgeqlf_f77 F77_FUNC(dgeqlf, DGEQLF)
#define dgeqp3_f77 F77_FUNC(dgeqp3, DGEQP3)
#define dgeqpf_f77 F77_FUNC(dgeqpf, DGEQPF)
#define dgeqr2_f77 F77_FUNC(dgeqr2, DGEQR2)
#define dgeqrf_f77 F77_FUNC(dgeqrf, DGEQRF)
#define dgerfs_f77 F77_FUNC(dgerfs, DGERFS)
#define dgerq2_f77 F77_FUNC(dgerq2, DGERQ2)
#define dgerqf_f77 F77_FUNC(dgerqf, DGERQF)
#define dgesc2_f77 F77_FUNC(dgesc2, DGESC2)
#define dgesdd_f77 F77_FUNC(dgesdd, DGESDD)
#define dgesv_f77 F77_FUNC(dgesv, DGESV)
#define dgesvd_f77 F77_FUNC(dgesvd, DGESVD)
#define dgesvx_f77 F77_FUNC(dgesvx, DGESVX)
#define dgetc2_f77 F77_FUNC(dgetc2, DGETC2)
#define dgetf2_f77 F77_FUNC(dgetf2, DGETF2)
#define dgetrf_f77 F77_FUNC(dgetrf, DGETRF)
#define dgetri_f77 F77_FUNC(dgetri, DGETRI)
#define dgetrs_f77 F77_FUNC(dgetrs, DGETRS)
#define dggbak_f77 F77_FUNC(dggbak, DGGBAK)
#define dggbal_f77 F77_FUNC(dggbal, DGGBAL)
#define dgges_f77 F77_FUNC(dgges, DGGES)
#define dggesx_f77 F77_FUNC(dggesx, DGGESX)
#define dggev_f77 F77_FUNC(dggev, DGGEV)
#define dggevx_f77 F77_FUNC(dggevx, DGGEVX)
#define dggglm_f77 F77_FUNC(dggglm, DGGGLM)
#define dgghrd_f77 F77_FUNC(dgghrd, DGGHRD)
#define dgglse_f77 F77_FUNC(dgglse, DGGLSE)
#define dggqrf_f77 F77_FUNC(dggqrf, DGGQRF)
#define dggrqf_f77 F77_FUNC(dggrqf, DGGRQF)
#define dggsvd_f77 F77_FUNC(dggsvd, DGGSVD)
#define dggsvp_f77 F77_FUNC(dggsvp, DGGSVP)
#define dgtcon_f77 F77_FUNC(dgtcon, DGTCON)
#define dgtrfs_f77 F77_FUNC(dgtrfs, DGTRFS)
#define dgtsv_f77 F77_FUNC(dgtsv, DGTSV)
#define dgtsvx_f77 F77_FUNC(dgtsvx, DGTSVX)
#define dgttrf_f77 F77_FUNC(dgttrf, DGTTRF)
#define dgttrs_f77 F77_FUNC(dgttrs, DGTTRS)
#define dgtts2_f77 F77_FUNC(dgtts2, DGTTS2)
#define dhgeqz_f77 F77_FUNC(dhgeqz, DHGEQZ)
#define dhsein_f77 F77_FUNC(dhsein, DHSEIN)
#define dhseqr_f77 F77_FUNC(dhseqr, DHSEQR)
#define disnan_f77 F77_FUNC(disnan, DISNAN)
#define dlabad_f77 F77_FUNC(dlabad, DLABAD)
#define dlabrd_f77 F77_FUNC(dlabrd, DLABRD)
#define dlacn2_f77 F77_FUNC(dlacn2, DLACN2)
#define dlacon_f77 F77_FUNC(dlacon, DLACON)
#define dlacpy_f77 F77_FUNC(dlacpy, DLACPY)
#define dladiv_f77 F77_FUNC(dladiv, DLADIV)
#define dlae2_f77 F77_FUNC(dlae2, DLAE2)
#define dlaebz_f77 F77_FUNC(dlaebz, DLAEBZ)
#define dlaed0_f77 F77_FUNC(dlaed0, DLAED0)
#define dlaed1_f77 F77_FUNC(dlaed1, DLAED1)
#define dlaed2_f77 F77_FUNC(dlaed2, DLAED2)
#define dlaed3_f77 F77_FUNC(dlaed3, DLAED3)
#define dlaed4_f77 F77_FUNC(dlaed4, DLAED4)
#define dlaed5_f77 F77_FUNC(dlaed5, DLAED5)
#define dlaed6_f77 F77_FUNC(dlaed6, DLAED6)
#define dlaed7_f77 F77_FUNC(dlaed7, DLAED7)
#define dlaed8_f77 F77_FUNC(dlaed8, DLAED8)
#define dlaed9_f77 F77_FUNC(dlaed9, DLAED9)
#define dlaeda_f77 F77_FUNC(dlaeda, DLAEDA)
#define dlaein_f77 F77_FUNC(dlaein, DLAEIN)
#define dlaev2_f77 F77_FUNC(dlaev2, DLAEV2)
#define dlaexc_f77 F77_FUNC(dlaexc, DLAEXC)
#define dlag2_f77 F77_FUNC(dlag2, DLAG2)
#define dlag2s_f77 F77_FUNC(dlag2s, DLAG2S)
#define dlags2_f77 F77_FUNC(dlags2, DLAGS2)
#define dlagtf_f77 F77_FUNC(dlagtf, DLAGTF)
#define dlagtm_f77 F77_FUNC(dlagtm, DLAGTM)
#define dlagts_f77 F77_FUNC(dlagts, DLAGTS)
#define dlagv2_f77 F77_FUNC(dlagv2, DLAGV2)
#define dlahqr_f77 F77_FUNC(dlahqr, DLAHQR)
#define dlahr2_f77 F77_FUNC(dlahr2, DLAHR2)
#define dlahrd_f77 F77_FUNC(dlahrd, DLAHRD)
#define dlaic1_f77 F77_FUNC(dlaic1, DLAIC1)
#define dlaisnan_f77 F77_FUNC(dlaisnan, DLAISNAN)
#define dlaln2_f77 F77_FUNC(dlaln2, DLALN2)
#define dlals0_f77 F77_FUNC(dlals0, DLALS0)
#define dlalsa_f77 F77_FUNC(dlalsa, DLALSA)
#define dlalsd_f77 F77_FUNC(dlalsd, DLALSD)
#define dlamch_f77 F77_FUNC(dlamch, DLAMCH)
#define dlamrg_f77 F77_FUNC(dlamrg, DLAMRG)
#define dlaneg_f77 F77_FUNC(dlaneg, DLANEG)
#define dlangb_f77 F77_FUNC(dlangb, DLANGB)
#define dlange_f77 F77_FUNC(dlange, DLANGE)
#define dlangt_f77 F77_FUNC(dlangt, DLANGT)
#define dlanhs_f77 F77_FUNC(dlanhs, DLANHS)
#define dlansb_f77 F77_FUNC(dlansb, DLANSB)
#define dlansp_f77 F77_FUNC(dlansp, DLANSP)
#define dlanst_f77 F77_FUNC(dlanst, DLANST)
#define dlansy_f77 F77_FUNC(dlansy, DLANSY)
#define dlantb_f77 F77_FUNC(dlantb, DLANTB)
#define dlantp_f77 F77_FUNC(dlantp, DLANTP)
#define dlantr_f77 F77_FUNC(dlantr, DLANTR)
#define dlanv2_f77 F77_FUNC(dlanv2, DLANV2)
#define dlapll_f77 F77_FUNC(dlapll, DLAPLL)
#define dlapmt_f77 F77_FUNC(dlapmt, DLAPMT)
#define dlapy2_f77 F77_FUNC(dlapy2, DLAPY2)
#define dlapy3_f77 F77_FUNC(dlapy3, DLAPY3)
#define dlaqgb_f77 F77_FUNC(dlaqgb, DLAQGB)
#define dlaqge_f77 F77_FUNC(dlaqge, DLAQGE)
#define dlaqp2_f77 F77_FUNC(dlaqp2, DLAQP2)
#define dlaqps_f77 F77_FUNC(dlaqps, DLAQPS)
#define dlaqr0_f77 F77_FUNC(dlaqr0, DLAQR0)
#define dlaqr1_f77 F77_FUNC(dlaqr1, DLAQR1)
#define dlaqr2_f77 F77_FUNC(dlaqr2, DLAQR2)
#define dlaqr3_f77 F77_FUNC(dlaqr3, DLAQR3)
#define dlaqr4_f77 F77_FUNC(dlaqr4, DLAQR4)
#define dlaqr5_f77 F77_FUNC(dlaqr5, DLAQR5)
#define dlaqsb_f77 F77_FUNC(dlaqsb, DLAQSB)
#define dlaqsp_f77 F77_FUNC(dlaqsp, DLAQSP)
#define dlaqsy_f77 F77_FUNC(dlaqsy, DLAQSY)
#define dlaqtr_f77 F77_FUNC(dlaqtr, DLAQTR)
#define dlar1v_f77 F77_FUNC(dlar1v, DLAR1V)
#define dlar2v_f77 F77_FUNC(dlar2v, DLAR2V)
#define dlarf_f77 F77_FUNC(dlarf, DLARF)
#define dlarfb_f77 F77_FUNC(dlarfb, DLARFB)
#define dlarfg_f77 F77_FUNC(dlarfg, DLARFG)
#define dlarft_f77 F77_FUNC(dlarft, DLARFT)
#define dlarfx_f77 F77_FUNC(dlarfx, DLARFX)
#define dlargv_f77 F77_FUNC(dlargv, DLARGV)
#define dlarnv_f77 F77_FUNC(dlarnv, DLARNV)
#define dlarra_f77 F77_FUNC(dlarra, DLARRA)
#define dlarrb_f77 F77_FUNC(dlarrb, DLARRB)
#define dlarrc_f77 F77_FUNC(dlarrc, DLARRC)
#define dlarrd_f77 F77_FUNC(dlarrd, DLARRD)
#define dlarre_f77 F77_FUNC(dlarre, DLARRE)
#define dlarrf_f77 F77_FUNC(dlarrf, DLARRF)
#define dlarrj_f77 F77_FUNC(dlarrj, DLARRJ)
#define dlarrk_f77 F77_FUNC(dlarrk, DLARRK)
#define dlarrr_f77 F77_FUNC(dlarrr, DLARRR)
#define dlarrv_f77 F77_FUNC(dlarrv, DLARRV)
#define dlartg_f77 F77_FUNC(dlartg, DLARTG)
#define dlartv_f77 F77_FUNC(dlartv, DLARTV)
#define dlaruv_f77 F77_FUNC(dlaruv, DLARUV)
#define dlarz_f77 F77_FUNC(dlarz, DLARZ)
#define dlarzb_f77 F77_FUNC(dlarzb, DLARZB)
#define dlarzt_f77 F77_FUNC(dlarzt, DLARZT)
#define dlas2_f77 F77_FUNC(dlas2, DLAS2)
#define dlascl_f77 F77_FUNC(dlascl, DLASCL)
#define dlasd0_f77 F77_FUNC(dlasd0, DLASD0)
#define dlasd1_f77 F77_FUNC(dlasd1, DLASD1)
#define dlasd2_f77 F77_FUNC(dlasd2, DLASD2)
#define dlasd3_f77 F77_FUNC(dlasd3, DLASD3)
#define dlasd4_f77 F77_FUNC(dlasd4, DLASD4)
#define dlasd5_f77 F77_FUNC(dlasd5, DLASD5)
#define dlasd6_f77 F77_FUNC(dlasd6, DLASD6)
#define dlasd7_f77 F77_FUNC(dlasd7, DLASD7)
#define dlasd8_f77 F77_FUNC(dlasd8, DLASD8)
#define dlasda_f77 F77_FUNC(dlasda, DLASDA)
#define dlasdq_f77 F77_FUNC(dlasdq, DLASDQ)
#define dlasdt_f77 F77_FUNC(dlasdt, DLASDT)
#define dlaset_f77 F77_FUNC(dlaset, DLASET)
#define dlasq1_f77 F77_FUNC(dlasq1, DLASQ1)
#define dlasq2_f77 F77_FUNC(dlasq2, DLASQ2)
#define dlasq3_f77 F77_FUNC(dlasq3, DLASQ3)
#define dlasq4_f77 F77_FUNC(dlasq4, DLASQ4)
#define dlasq5_f77 F77_FUNC(dlasq5, DLASQ5)
#define dlasq6_f77 F77_FUNC(dlasq6, DLASQ6)
#define dlasr_f77 F77_FUNC(dlasr, DLASR)
#define dlasrt_f77 F77_FUNC(dlasrt, DLASRT)
#define dlassq_f77 F77_FUNC(dlassq, DLASSQ)
#define dlasv2_f77 F77_FUNC(dlasv2, DLASV2)
#define dlaswp_f77 F77_FUNC(dlaswp, DLASWP)
#define dlasy2_f77 F77_FUNC(dlasy2, DLASY2)
#define dlasyf_f77 F77_FUNC(dlasyf, DLASYF)
#define dlatbs_f77 F77_FUNC(dlatbs, DLATBS)
#define dlatdf_f77 F77_FUNC(dlatdf, DLATDF)
#define dlatps_f77 F77_FUNC(dlatps, DLATPS)
#define dlatrd_f77 F77_FUNC(dlatrd, DLATRD)
#define dlatrs_f77 F77_FUNC(dlatrs, DLATRS)
#define dlatrz_f77 F77_FUNC(dlatrz, DLATRZ)
#define dlatzm_f77 F77_FUNC(dlatzm, DLATZM)
#define dlauu2_f77 F77_FUNC(dlauu2, DLAUU2)
#define dlauum_f77 F77_FUNC(dlauum, DLAUUM)
#define dlazq3_f77 F77_FUNC(dlazq3, DLAZQ3)
#define dlazq4_f77 F77_FUNC(dlazq4, DLAZQ4)
#define dopgtr_f77 F77_FUNC(dopgtr, DOPGTR)
#define dopmtr_f77 F77_FUNC(dopmtr, DOPMTR)
#define dorg2l_f77 F77_FUNC(dorg2l, DORG2L)
#define dorg2r_f77 F77_FUNC(dorg2r, DORG2R)
#define dorgbr_f77 F77_FUNC(dorgbr, DORGBR)
#define dorghr_f77 F77_FUNC(dorghr, DORGHR)
#define dorgl2_f77 F77_FUNC(dorgl2, DORGL2)
#define dorglq_f77 F77_FUNC(dorglq, DORGLQ)
#define dorgql_f77 F77_FUNC(dorgql, DORGQL)
#define dorgqr_f77 F77_FUNC(dorgqr, DORGQR)
#define dorgr2_f77 F77_FUNC(dorgr2, DORGR2)
#define dorgrq_f77 F77_FUNC(dorgrq, DORGRQ)
#define dorgtr_f77 F77_FUNC(dorgtr, DORGTR)
#define dorm2l_f77 F77_FUNC(dorm2l, DORM2L)
#define dorm2r_f77 F77_FUNC(dorm2r, DORM2R)
#define dormbr_f77 F77_FUNC(dormbr, DORMBR)
#define dormhr_f77 F77_FUNC(dormhr, DORMHR)
#define dorml2_f77 F77_FUNC(dorml2, DORML2)
#define dormlq_f77 F77_FUNC(dormlq, DORMLQ)
#define dormql_f77 F77_FUNC(dormql, DORMQL)
#define dormqr_f77 F77_FUNC(dormqr, DORMQR)
#define dormr2_f77 F77_FUNC(dormr2, DORMR2)
#define dormr3_f77 F77_FUNC(dormr3, DORMR3)
#define dormrq_f77 F77_FUNC(dormrq, DORMRQ)
#define dormrz_f77 F77_FUNC(dormrz, DORMRZ)
#define dormtr_f77 F77_FUNC(dormtr, DORMTR)
#define dpbcon_f77 F77_FUNC(dpbcon, DPBCON)
#define dpbequ_f77 F77_FUNC(dpbequ, DPBEQU)
#define dpbrfs_f77 F77_FUNC(dpbrfs, DPBRFS)
#define dpbstf_f77 F77_FUNC(dpbstf, DPBSTF)
#define dpbsv_f77 F77_FUNC(dpbsv, DPBSV)
#define dpbsvx_f77 F77_FUNC(dpbsvx, DPBSVX)
#define dpbtf2_f77 F77_FUNC(dpbtf2, DPBTF2)
#define dpbtrf_f77 F77_FUNC(dpbtrf, DPBTRF)
#define dpbtrs_f77 F77_FUNC(dpbtrs, DPBTRS)
#define dpocon_f77 F77_FUNC(dpocon, DPOCON)
#define dpoequ_f77 F77_FUNC(dpoequ, DPOEQU)
#define dporfs_f77 F77_FUNC(dporfs, DPORFS)
#define dposv_f77 F77_FUNC(dposv, DPOSV)
#define dposvx_f77 F77_FUNC(dposvx, DPOSVX)
#define dpotf2_f77 F77_FUNC(dpotf2, DPOTF2)
#define dpotrf_f77 F77_FUNC(dpotrf, DPOTRF)
#define dpotri_f77 F77_FUNC(dpotri, DPOTRI)
#define dpotrs_f77 F77_FUNC(dpotrs, DPOTRS)
#define dppcon_f77 F77_FUNC(dppcon, DPPCON)
#define dppequ_f77 F77_FUNC(dppequ, DPPEQU)
#define dpprfs_f77 F77_FUNC(dpprfs, DPPRFS)
#define dppsv_f77 F77_FUNC(dppsv, DPPSV)
#define dppsvx_f77 F77_FUNC(dppsvx, DPPSVX)
#define dpptrf_f77 F77_FUNC(dpptrf, DPPTRF)
#define dpptri_f77 F77_FUNC(dpptri, DPPTRI)
#define dpptrs_f77 F77_FUNC(dpptrs, DPPTRS)
#define dptcon_f77 F77_FUNC(dptcon, DPTCON)
#define dpteqr_f77 F77_FUNC(dpteqr, DPTEQR)
#define dptrfs_f77 F77_FUNC(dptrfs, DPTRFS)
#define dptsv_f77 F77_FUNC(dptsv, DPTSV)
#define dptsvx_f77 F77_FUNC(dptsvx, DPTSVX)
#define dpttrf_f77 F77_FUNC(dpttrf, DPTTRF)
#define dpttrs_f77 F77_FUNC(dpttrs, DPTTRS)
#define dptts2_f77 F77_FUNC(dptts2, DPTTS2)
#define drscl_f77 F77_FUNC(drscl, DRSCL)
#define dsbev_f77 F77_FUNC(dsbev, DSBEV)
#define dsbevd_f77 F77_FUNC(dsbevd, DSBEVD)
#define dsbevx_f77 F77_FUNC(dsbevx, DSBEVX)
#define dsbgst_f77 F77_FUNC(dsbgst, DSBGST)
#define dsbgv_f77 F77_FUNC(dsbgv, DSBGV)
#define dsbgvd_f77 F77_FUNC(dsbgvd, DSBGVD)
#define dsbgvx_f77 F77_FUNC(dsbgvx, DSBGVX)
#define dsbtrd_f77 F77_FUNC(dsbtrd, DSBTRD)
#define dsgesv_f77 F77_FUNC(dsgesv, DSGESV)
#define dspcon_f77 F77_FUNC(dspcon, DSPCON)
#define dspev_f77 F77_FUNC(dspev, DSPEV)
#define dspevd_f77 F77_FUNC(dspevd, DSPEVD)
#define dspevx_f77 F77_FUNC(dspevx, DSPEVX)
#define dspgst_f77 F77_FUNC(dspgst, DSPGST)
#define dspgv_f77 F77_FUNC(dspgv, DSPGV)
#define dspgvd_f77 F77_FUNC(dspgvd, DSPGVD)
#define dspgvx_f77 F77_FUNC(dspgvx, DSPGVX)
#define dsprfs_f77 F77_FUNC(dsprfs, DSPRFS)
#define dspsv_f77 F77_FUNC(dspsv, DSPSV)
#define dspsvx_f77 F77_FUNC(dspsvx, DSPSVX)
#define dsptrd_f77 F77_FUNC(dsptrd, DSPTRD)
#define dsptrf_f77 F77_FUNC(dsptrf, DSPTRF)
#define dsptri_f77 F77_FUNC(dsptri, DSPTRI)
#define dsptrs_f77 F77_FUNC(dsptrs, DSPTRS)
#define dstebz_f77 F77_FUNC(dstebz, DSTEBZ)
#define dstedc_f77 F77_FUNC(dstedc, DSTEDC)
#define dstegr_f77 F77_FUNC(dstegr, DSTEGR)
#define dstein_f77 F77_FUNC(dstein, DSTEIN)
#define dstemr_f77 F77_FUNC(dstemr, DSTEMR)
#define dsteqr_f77 F77_FUNC(dsteqr, DSTEQR)
#define dsterf_f77 F77_FUNC(dsterf, DSTERF)
#define dstev_f77 F77_FUNC(dstev, DSTEV)
#define dstevd_f77 F77_FUNC(dstevd, DSTEVD)
#define dstevr_f77 F77_FUNC(dstevr, DSTEVR)
#define dstevx_f77 F77_FUNC(dstevx, DSTEVX)
#define dsycon_f77 F77_FUNC(dsycon, DSYCON)
#define dsyev_f77 F77_FUNC(dsyev, DSYEV)
#define dsyevd_f77 F77_FUNC(dsyevd, DSYEVD)
#define dsyevr_f77 F77_FUNC(dsyevr, DSYEVR)
#define dsyevx_f77 F77_FUNC(dsyevx, DSYEVX)
#define dsygs2_f77 F77_FUNC(dsygs2, DSYGS2)
#define dsygst_f77 F77_FUNC(dsygst, DSYGST)
#define dsygv_f77 F77_FUNC(dsygv, DSYGV)
#define dsygvd_f77 F77_FUNC(dsygvd, DSYGVD)
#define dsygvx_f77 F77_FUNC(dsygvx, DSYGVX)
#define dsyrfs_f77 F77_FUNC(dsyrfs, DSYRFS)
#define dsysv_f77 F77_FUNC(dsysv, DSYSV)
#define dsysvx_f77 F77_FUNC(dsysvx, DSYSVX)
#define dsytd2_f77 F77_FUNC(dsytd2, DSYTD2)
#define dsytf2_f77 F77_FUNC(dsytf2, DSYTF2)
#define dsytrd_f77 F77_FUNC(dsytrd, DSYTRD)
#define dsytrf_f77 F77_FUNC(dsytrf, DSYTRF)
#define dsytri_f77 F77_FUNC(dsytri, DSYTRI)
#define dsytrs_f77 F77_FUNC(dsytrs, DSYTRS)
#define dtbcon_f77 F77_FUNC(dtbcon, DTBCON)
#define dtbrfs_f77 F77_FUNC(dtbrfs, DTBRFS)
#define dtbtrs_f77 F77_FUNC(dtbtrs, DTBTRS)
#define dtgevc_f77 F77_FUNC(dtgevc, DTGEVC)
#define dtgex2_f77 F77_FUNC(dtgex2, DTGEX2)
#define dtgexc_f77 F77_FUNC(dtgexc, DTGEXC)
#define dtgsen_f77 F77_FUNC(dtgsen, DTGSEN)
#define dtgsja_f77 F77_FUNC(dtgsja, DTGSJA)
#define dtgsna_f77 F77_FUNC(dtgsna, DTGSNA)
#define dtgsy2_f77 F77_FUNC(dtgsy2, DTGSY2)
#define dtgsyl_f77 F77_FUNC(dtgsyl, DTGSYL)
#define dtpcon_f77 F77_FUNC(dtpcon, DTPCON)
#define dtprfs_f77 F77_FUNC(dtprfs, DTPRFS)
#define dtptri_f77 F77_FUNC(dtptri, DTPTRI)
#define dtptrs_f77 F77_FUNC(dtptrs, DTPTRS)
#define dtrcon_f77 F77_FUNC(dtrcon, DTRCON)
#define dtrevc_f77 F77_FUNC(dtrevc, DTREVC)
#define dtrexc_f77 F77_FUNC(dtrexc, DTREXC)
#define dtrrfs_f77 F77_FUNC(dtrrfs, DTRRFS)
#define dtrsen_f77 F77_FUNC(dtrsen, DTRSEN)
#define dtrsna_f77 F77_FUNC(dtrsna, DTRSNA)
#define dtrsyl_f77 F77_FUNC(dtrsyl, DTRSYL)
#define dtrti2_f77 F77_FUNC(dtrti2, DTRTI2)
#define dtrtri_f77 F77_FUNC(dtrtri, DTRTRI)
#define dtrtrs_f77 F77_FUNC(dtrtrs, DTRTRS)
#define dtzrqf_f77 F77_FUNC(dtzrqf, DTZRQF)
#define dtzrzf_f77 F77_FUNC(dtzrzf, DTZRZF)
#define dzsum1_f77 F77_FUNC(dzsum1, DZSUM1)
#define icmax1_f77 F77_FUNC(icmax1, ICMAX1)
#define ieeeck_f77 F77_FUNC(ieeeck, IEEECK)
#define ilaenv_f77 F77_FUNC(ilaenv, ILAENV)
#define ilaver_f77 F77_FUNC(ilaver, ILAVER)
#define iparmq_f77 F77_FUNC(iparmq, IPARMQ)
#define izmax1_f77 F77_FUNC(izmax1, IZMAX1)
#define lsamen_f77 F77_FUNC(lsamen, LSAMEN)
#define sbdsdc_f77 F77_FUNC(sbdsdc, SBDSDC)
#define sbdsqr_f77 F77_FUNC(sbdsqr, SBDSQR)
#define scsum1_f77 F77_FUNC(scsum1, SCSUM1)
#define sdisna_f77 F77_FUNC(sdisna, SDISNA)
#define sgbbrd_f77 F77_FUNC(sgbbrd, SGBBRD)
#define sgbcon_f77 F77_FUNC(sgbcon, SGBCON)
#define sgbequ_f77 F77_FUNC(sgbequ, SGBEQU)
#define sgbrfs_f77 F77_FUNC(sgbrfs, SGBRFS)
#define sgbsv_f77 F77_FUNC(sgbsv, SGBSV)
#define sgbsvx_f77 F77_FUNC(sgbsvx, SGBSVX)
#define sgbtf2_f77 F77_FUNC(sgbtf2, SGBTF2)
#define sgbtrf_f77 F77_FUNC(sgbtrf, SGBTRF)
#define sgbtrs_f77 F77_FUNC(sgbtrs, SGBTRS)
#define sgebak_f77 F77_FUNC(sgebak, SGEBAK)
#define sgebal_f77 F77_FUNC(sgebal, SGEBAL)
#define sgebd2_f77 F77_FUNC(sgebd2, SGEBD2)
#define sgebrd_f77 F77_FUNC(sgebrd, SGEBRD)
#define sgecon_f77 F77_FUNC(sgecon, SGECON)
#define sgeequ_f77 F77_FUNC(sgeequ, SGEEQU)
#define sgees_f77 F77_FUNC(sgees, SGEES)
#define sgeesx_f77 F77_FUNC(sgeesx, SGEESX)
#define sgeev_f77 F77_FUNC(sgeev, SGEEV)
#define sgeevx_f77 F77_FUNC(sgeevx, SGEEVX)
#define sgegs_f77 F77_FUNC(sgegs, SGEGS)
#define sgegv_f77 F77_FUNC(sgegv, SGEGV)
#define sgehd2_f77 F77_FUNC(sgehd2, SGEHD2)
#define sgehrd_f77 F77_FUNC(sgehrd, SGEHRD)
#define sgelq2_f77 F77_FUNC(sgelq2, SGELQ2)
#define sgelqf_f77 F77_FUNC(sgelqf, SGELQF)
#define sgels_f77 F77_FUNC(sgels, SGELS)
#define sgelsd_f77 F77_FUNC(sgelsd, SGELSD)
#define sgelss_f77 F77_FUNC(sgelss, SGELSS)
#define sgelsx_f77 F77_FUNC(sgelsx, SGELSX)
#define sgelsy_f77 F77_FUNC(sgelsy, SGELSY)
#define sgeql2_f77 F77_FUNC(sgeql2, SGEQL2)
#define sgeqlf_f77 F77_FUNC(sgeqlf, SGEQLF)
#define sgeqp3_f77 F77_FUNC(sgeqp3, SGEQP3)
#define sgeqpf_f77 F77_FUNC(sgeqpf, SGEQPF)
#define sgeqr2_f77 F77_FUNC(sgeqr2, SGEQR2)
#define sgeqrf_f77 F77_FUNC(sgeqrf, SGEQRF)
#define sgerfs_f77 F77_FUNC(sgerfs, SGERFS)
#define sgerq2_f77 F77_FUNC(sgerq2, SGERQ2)
#define sgerqf_f77 F77_FUNC(sgerqf, SGERQF)
#define sgesc2_f77 F77_FUNC(sgesc2, SGESC2)
#define sgesdd_f77 F77_FUNC(sgesdd, SGESDD)
#define sgesv_f77 F77_FUNC(sgesv, SGESV)
#define sgesvd_f77 F77_FUNC(sgesvd, SGESVD)
#define sgesvx_f77 F77_FUNC(sgesvx, SGESVX)
#define sgetc2_f77 F77_FUNC(sgetc2, SGETC2)
#define sgetf2_f77 F77_FUNC(sgetf2, SGETF2)
#define sgetrf_f77 F77_FUNC(sgetrf, SGETRF)
#define sgetri_f77 F77_FUNC(sgetri, SGETRI)
#define sgetrs_f77 F77_FUNC(sgetrs, SGETRS)
#define sggbak_f77 F77_FUNC(sggbak, SGGBAK)
#define sggbal_f77 F77_FUNC(sggbal, SGGBAL)
#define sgges_f77 F77_FUNC(sgges, SGGES)
#define sggesx_f77 F77_FUNC(sggesx, SGGESX)
#define sggev_f77 F77_FUNC(sggev, SGGEV)
#define sggevx_f77 F77_FUNC(sggevx, SGGEVX)
#define sggglm_f77 F77_FUNC(sggglm, SGGGLM)
#define sgghrd_f77 F77_FUNC(sgghrd, SGGHRD)
#define sgglse_f77 F77_FUNC(sgglse, SGGLSE)
#define sggqrf_f77 F77_FUNC(sggqrf, SGGQRF)
#define sggrqf_f77 F77_FUNC(sggrqf, SGGRQF)
#define sggsvd_f77 F77_FUNC(sggsvd, SGGSVD)
#define sggsvp_f77 F77_FUNC(sggsvp, SGGSVP)
#define sgtcon_f77 F77_FUNC(sgtcon, SGTCON)
#define sgtrfs_f77 F77_FUNC(sgtrfs, SGTRFS)
#define sgtsv_f77 F77_FUNC(sgtsv, SGTSV)
#define sgtsvx_f77 F77_FUNC(sgtsvx, SGTSVX)
#define sgttrf_f77 F77_FUNC(sgttrf, SGTTRF)
#define sgttrs_f77 F77_FUNC(sgttrs, SGTTRS)
#define sgtts2_f77 F77_FUNC(sgtts2, SGTTS2)
#define shgeqz_f77 F77_FUNC(shgeqz, SHGEQZ)
#define shsein_f77 F77_FUNC(shsein, SHSEIN)
#define shseqr_f77 F77_FUNC(shseqr, SHSEQR)
#define sisnan_f77 F77_FUNC(sisnan, SISNAN)
#define slabad_f77 F77_FUNC(slabad, SLABAD)
#define slabrd_f77 F77_FUNC(slabrd, SLABRD)
#define slacn2_f77 F77_FUNC(slacn2, SLACN2)
#define slacon_f77 F77_FUNC(slacon, SLACON)
#define slacpy_f77 F77_FUNC(slacpy, SLACPY)
#define sladiv_f77 F77_FUNC(sladiv, SLADIV)
#define slae2_f77 F77_FUNC(slae2, SLAE2)
#define slaebz_f77 F77_FUNC(slaebz, SLAEBZ)
#define slaed0_f77 F77_FUNC(slaed0, SLAED0)
#define slaed1_f77 F77_FUNC(slaed1, SLAED1)
#define slaed2_f77 F77_FUNC(slaed2, SLAED2)
#define slaed3_f77 F77_FUNC(slaed3, SLAED3)
#define slaed4_f77 F77_FUNC(slaed4, SLAED4)
#define slaed5_f77 F77_FUNC(slaed5, SLAED5)
#define slaed6_f77 F77_FUNC(slaed6, SLAED6)
#define slaed7_f77 F77_FUNC(slaed7, SLAED7)
#define slaed8_f77 F77_FUNC(slaed8, SLAED8)
#define slaed9_f77 F77_FUNC(slaed9, SLAED9)
#define slaeda_f77 F77_FUNC(slaeda, SLAEDA)
#define slaein_f77 F77_FUNC(slaein, SLAEIN)
#define slaev2_f77 F77_FUNC(slaev2, SLAEV2)
#define slaexc_f77 F77_FUNC(slaexc, SLAEXC)
#define slag2_f77 F77_FUNC(slag2, SLAG2)
#define slag2d_f77 F77_FUNC(slag2d, SLAG2D)
#define slags2_f77 F77_FUNC(slags2, SLAGS2)
#define slagtf_f77 F77_FUNC(slagtf, SLAGTF)
#define slagtm_f77 F77_FUNC(slagtm, SLAGTM)
#define slagts_f77 F77_FUNC(slagts, SLAGTS)
#define slagv2_f77 F77_FUNC(slagv2, SLAGV2)
#define slahqr_f77 F77_FUNC(slahqr, SLAHQR)
#define slahr2_f77 F77_FUNC(slahr2, SLAHR2)
#define slahrd_f77 F77_FUNC(slahrd, SLAHRD)
#define slaic1_f77 F77_FUNC(slaic1, SLAIC1)
#define slaisnan_f77 F77_FUNC(slaisnan, SLAISNAN)
#define slaln2_f77 F77_FUNC(slaln2, SLALN2)
#define slals0_f77 F77_FUNC(slals0, SLALS0)
#define slalsa_f77 F77_FUNC(slalsa, SLALSA)
#define slalsd_f77 F77_FUNC(slalsd, SLALSD)
#define slamrg_f77 F77_FUNC(slamrg, SLAMRG)
#define slaneg_f77 F77_FUNC(slaneg, SLANEG)
#define slangb_f77 F77_FUNC(slangb, SLANGB)
#define slange_f77 F77_FUNC(slange, SLANGE)
#define slangt_f77 F77_FUNC(slangt, SLANGT)
#define slanhs_f77 F77_FUNC(slanhs, SLANHS)
#define slansb_f77 F77_FUNC(slansb, SLANSB)
#define slansp_f77 F77_FUNC(slansp, SLANSP)
#define slanst_f77 F77_FUNC(slanst, SLANST)
#define slansy_f77 F77_FUNC(slansy, SLANSY)
#define slantb_f77 F77_FUNC(slantb, SLANTB)
#define slantp_f77 F77_FUNC(slantp, SLANTP)
#define slantr_f77 F77_FUNC(slantr, SLANTR)
#define slanv2_f77 F77_FUNC(slanv2, SLANV2)
#define slapll_f77 F77_FUNC(slapll, SLAPLL)
#define slapmt_f77 F77_FUNC(slapmt, SLAPMT)
#define slapy2_f77 F77_FUNC(slapy2, SLAPY2)
#define slapy3_f77 F77_FUNC(slapy3, SLAPY3)
#define slaqgb_f77 F77_FUNC(slaqgb, SLAQGB)
#define slaqge_f77 F77_FUNC(slaqge, SLAQGE)
#define slaqp2_f77 F77_FUNC(slaqp2, SLAQP2)
#define slaqps_f77 F77_FUNC(slaqps, SLAQPS)
#define slaqr0_f77 F77_FUNC(slaqr0, SLAQR0)
#define slaqr1_f77 F77_FUNC(slaqr1, SLAQR1)
#define slaqr2_f77 F77_FUNC(slaqr2, SLAQR2)
#define slaqr3_f77 F77_FUNC(slaqr3, SLAQR3)
#define slaqr4_f77 F77_FUNC(slaqr4, SLAQR4)
#define slaqr5_f77 F77_FUNC(slaqr5, SLAQR5)
#define slaqsb_f77 F77_FUNC(slaqsb, SLAQSB)
#define slaqsp_f77 F77_FUNC(slaqsp, SLAQSP)
#define slaqsy_f77 F77_FUNC(slaqsy, SLAQSY)
#define slaqtr_f77 F77_FUNC(slaqtr, SLAQTR)
#define slar1v_f77 F77_FUNC(slar1v, SLAR1V)
#define slar2v_f77 F77_FUNC(slar2v, SLAR2V)
#define slarf_f77 F77_FUNC(slarf, SLARF)
#define slarfb_f77 F77_FUNC(slarfb, SLARFB)
#define slarfg_f77 F77_FUNC(slarfg, SLARFG)
#define slarft_f77 F77_FUNC(slarft, SLARFT)
#define slarfx_f77 F77_FUNC(slarfx, SLARFX)
#define slargv_f77 F77_FUNC(slargv, SLARGV)
#define slarnv_f77 F77_FUNC(slarnv, SLARNV)
#define slarra_f77 F77_FUNC(slarra, SLARRA)
#define slarrb_f77 F77_FUNC(slarrb, SLARRB)
#define slarrc_f77 F77_FUNC(slarrc, SLARRC)
#define slarrd_f77 F77_FUNC(slarrd, SLARRD)
#define slarre_f77 F77_FUNC(slarre, SLARRE)
#define slarrf_f77 F77_FUNC(slarrf, SLARRF)
#define slarrj_f77 F77_FUNC(slarrj, SLARRJ)
#define slarrk_f77 F77_FUNC(slarrk, SLARRK)
#define slarrr_f77 F77_FUNC(slarrr, SLARRR)
#define slarrv_f77 F77_FUNC(slarrv, SLARRV)
#define slartg_f77 F77_FUNC(slartg, SLARTG)
#define slartv_f77 F77_FUNC(slartv, SLARTV)
#define slaruv_f77 F77_FUNC(slaruv, SLARUV)
#define slarz_f77 F77_FUNC(slarz, SLARZ)
#define slarzb_f77 F77_FUNC(slarzb, SLARZB)
#define slarzt_f77 F77_FUNC(slarzt, SLARZT)
#define slas2_f77 F77_FUNC(slas2, SLAS2)
#define slascl_f77 F77_FUNC(slascl, SLASCL)
#define slasd0_f77 F77_FUNC(slasd0, SLASD0)
#define slasd1_f77 F77_FUNC(slasd1, SLASD1)
#define slasd2_f77 F77_FUNC(slasd2, SLASD2)
#define slasd3_f77 F77_FUNC(slasd3, SLASD3)
#define slasd4_f77 F77_FUNC(slasd4, SLASD4)
#define slasd5_f77 F77_FUNC(slasd5, SLASD5)
#define slasd6_f77 F77_FUNC(slasd6, SLASD6)
#define slasd7_f77 F77_FUNC(slasd7, SLASD7)
#define slasd8_f77 F77_FUNC(slasd8, SLASD8)
#define slasda_f77 F77_FUNC(slasda, SLASDA)
#define slasdq_f77 F77_FUNC(slasdq, SLASDQ)
#define slasdt_f77 F77_FUNC(slasdt, SLASDT)
#define slaset_f77 F77_FUNC(slaset, SLASET)
#define slasq1_f77 F77_FUNC(slasq1, SLASQ1)
#define slasq2_f77 F77_FUNC(slasq2, SLASQ2)
#define slasq3_f77 F77_FUNC(slasq3, SLASQ3)
#define slasq4_f77 F77_FUNC(slasq4, SLASQ4)
#define slasq5_f77 F77_FUNC(slasq5, SLASQ5)
#define slasq6_f77 F77_FUNC(slasq6, SLASQ6)
#define slasr_f77 F77_FUNC(slasr, SLASR)
#define slasrt_f77 F77_FUNC(slasrt, SLASRT)
#define slassq_f77 F77_FUNC(slassq, SLASSQ)
#define slasv2_f77 F77_FUNC(slasv2, SLASV2)
#define slaswp_f77 F77_FUNC(slaswp, SLASWP)
#define slasy2_f77 F77_FUNC(slasy2, SLASY2)
#define slasyf_f77 F77_FUNC(slasyf, SLASYF)
#define slatbs_f77 F77_FUNC(slatbs, SLATBS)
#define slatdf_f77 F77_FUNC(slatdf, SLATDF)
#define slatps_f77 F77_FUNC(slatps, SLATPS)
#define slatrd_f77 F77_FUNC(slatrd, SLATRD)
#define slatrs_f77 F77_FUNC(slatrs, SLATRS)
#define slatrz_f77 F77_FUNC(slatrz, SLATRZ)
#define slatzm_f77 F77_FUNC(slatzm, SLATZM)
#define slauu2_f77 F77_FUNC(slauu2, SLAUU2)
#define slauum_f77 F77_FUNC(slauum, SLAUUM)
#define slazq3_f77 F77_FUNC(slazq3, SLAZQ3)
#define slazq4_f77 F77_FUNC(slazq4, SLAZQ4)
#define sopgtr_f77 F77_FUNC(sopgtr, SOPGTR)
#define sopmtr_f77 F77_FUNC(sopmtr, SOPMTR)
#define sorg2l_f77 F77_FUNC(sorg2l, SORG2L)
#define sorg2r_f77 F77_FUNC(sorg2r, SORG2R)
#define sorgbr_f77 F77_FUNC(sorgbr, SORGBR)
#define sorghr_f77 F77_FUNC(sorghr, SORGHR)
#define sorgl2_f77 F77_FUNC(sorgl2, SORGL2)
#define sorglq_f77 F77_FUNC(sorglq, SORGLQ)
#define sorgql_f77 F77_FUNC(sorgql, SORGQL)
#define sorgqr_f77 F77_FUNC(sorgqr, SORGQR)
#define sorgr2_f77 F77_FUNC(sorgr2, SORGR2)
#define sorgrq_f77 F77_FUNC(sorgrq, SORGRQ)
#define sorgtr_f77 F77_FUNC(sorgtr, SORGTR)
#define sorm2l_f77 F77_FUNC(sorm2l, SORM2L)
#define sorm2r_f77 F77_FUNC(sorm2r, SORM2R)
#define sormbr_f77 F77_FUNC(sormbr, SORMBR)
#define sormhr_f77 F77_FUNC(sormhr, SORMHR)
#define sorml2_f77 F77_FUNC(sorml2, SORML2)
#define sormlq_f77 F77_FUNC(sormlq, SORMLQ)
#define sormql_f77 F77_FUNC(sormql, SORMQL)
#define sormqr_f77 F77_FUNC(sormqr, SORMQR)
#define sormr2_f77 F77_FUNC(sormr2, SORMR2)
#define sormr3_f77 F77_FUNC(sormr3, SORMR3)
#define sormrq_f77 F77_FUNC(sormrq, SORMRQ)
#define sormrz_f77 F77_FUNC(sormrz, SORMRZ)
#define sormtr_f77 F77_FUNC(sormtr, SORMTR)
#define spbcon_f77 F77_FUNC(spbcon, SPBCON)
#define spbequ_f77 F77_FUNC(spbequ, SPBEQU)
#define spbrfs_f77 F77_FUNC(spbrfs, SPBRFS)
#define spbstf_f77 F77_FUNC(spbstf, SPBSTF)
#define spbsv_f77 F77_FUNC(spbsv, SPBSV)
#define spbsvx_f77 F77_FUNC(spbsvx, SPBSVX)
#define spbtf2_f77 F77_FUNC(spbtf2, SPBTF2)
#define spbtrf_f77 F77_FUNC(spbtrf, SPBTRF)
#define spbtrs_f77 F77_FUNC(spbtrs, SPBTRS)
#define spocon_f77 F77_FUNC(spocon, SPOCON)
#define spoequ_f77 F77_FUNC(spoequ, SPOEQU)
#define sporfs_f77 F77_FUNC(sporfs, SPORFS)
#define sposv_f77 F77_FUNC(sposv, SPOSV)
#define sposvx_f77 F77_FUNC(sposvx, SPOSVX)
#define spotf2_f77 F77_FUNC(spotf2, SPOTF2)
#define spotrf_f77 F77_FUNC(spotrf, SPOTRF)
#define spotri_f77 F77_FUNC(spotri, SPOTRI)
#define spotrs_f77 F77_FUNC(spotrs, SPOTRS)
#define sppcon_f77 F77_FUNC(sppcon, SPPCON)
#define sppequ_f77 F77_FUNC(sppequ, SPPEQU)
#define spprfs_f77 F77_FUNC(spprfs, SPPRFS)
#define sppsv_f77 F77_FUNC(sppsv, SPPSV)
#define sppsvx_f77 F77_FUNC(sppsvx, SPPSVX)
#define spptrf_f77 F77_FUNC(spptrf, SPPTRF)
#define spptri_f77 F77_FUNC(spptri, SPPTRI)
#define spptrs_f77 F77_FUNC(spptrs, SPPTRS)
#define sptcon_f77 F77_FUNC(sptcon, SPTCON)
#define spteqr_f77 F77_FUNC(spteqr, SPTEQR)
#define sptrfs_f77 F77_FUNC(sptrfs, SPTRFS)
#define sptsv_f77 F77_FUNC(sptsv, SPTSV)
#define sptsvx_f77 F77_FUNC(sptsvx, SPTSVX)
#define spttrf_f77 F77_FUNC(spttrf, SPTTRF)
#define spttrs_f77 F77_FUNC(spttrs, SPTTRS)
#define sptts2_f77 F77_FUNC(sptts2, SPTTS2)
#define srscl_f77 F77_FUNC(srscl, SRSCL)
#define ssbev_f77 F77_FUNC(ssbev, SSBEV)
#define ssbevd_f77 F77_FUNC(ssbevd, SSBEVD)
#define ssbevx_f77 F77_FUNC(ssbevx, SSBEVX)
#define ssbgst_f77 F77_FUNC(ssbgst, SSBGST)
#define ssbgv_f77 F77_FUNC(ssbgv, SSBGV)
#define ssbgvd_f77 F77_FUNC(ssbgvd, SSBGVD)
#define ssbgvx_f77 F77_FUNC(ssbgvx, SSBGVX)
#define ssbtrd_f77 F77_FUNC(ssbtrd, SSBTRD)
#define sspcon_f77 F77_FUNC(sspcon, SSPCON)
#define sspev_f77 F77_FUNC(sspev, SSPEV)
#define sspevd_f77 F77_FUNC(sspevd, SSPEVD)
#define sspevx_f77 F77_FUNC(sspevx, SSPEVX)
#define sspgst_f77 F77_FUNC(sspgst, SSPGST)
#define sspgv_f77 F77_FUNC(sspgv, SSPGV)
#define sspgvd_f77 F77_FUNC(sspgvd, SSPGVD)
#define sspgvx_f77 F77_FUNC(sspgvx, SSPGVX)
#define ssprfs_f77 F77_FUNC(ssprfs, SSPRFS)
#define sspsv_f77 F77_FUNC(sspsv, SSPSV)
#define sspsvx_f77 F77_FUNC(sspsvx, SSPSVX)
#define ssptrd_f77 F77_FUNC(ssptrd, SSPTRD)
#define ssptrf_f77 F77_FUNC(ssptrf, SSPTRF)
#define ssptri_f77 F77_FUNC(ssptri, SSPTRI)
#define ssptrs_f77 F77_FUNC(ssptrs, SSPTRS)
#define sstebz_f77 F77_FUNC(sstebz, SSTEBZ)
#define sstedc_f77 F77_FUNC(sstedc, SSTEDC)
#define sstegr_f77 F77_FUNC(sstegr, SSTEGR)
#define sstein_f77 F77_FUNC(sstein, SSTEIN)
#define sstemr_f77 F77_FUNC(sstemr, SSTEMR)
#define ssteqr_f77 F77_FUNC(ssteqr, SSTEQR)
#define ssterf_f77 F77_FUNC(ssterf, SSTERF)
#define sstev_f77 F77_FUNC(sstev, SSTEV)
#define sstevd_f77 F77_FUNC(sstevd, SSTEVD)
#define sstevr_f77 F77_FUNC(sstevr, SSTEVR)
#define sstevx_f77 F77_FUNC(sstevx, SSTEVX)
#define ssycon_f77 F77_FUNC(ssycon, SSYCON)
#define ssyev_f77 F77_FUNC(ssyev, SSYEV)
#define ssyevd_f77 F77_FUNC(ssyevd, SSYEVD)
#define ssyevr_f77 F77_FUNC(ssyevr, SSYEVR)
#define ssyevx_f77 F77_FUNC(ssyevx, SSYEVX)
#define ssygs2_f77 F77_FUNC(ssygs2, SSYGS2)
#define ssygst_f77 F77_FUNC(ssygst, SSYGST)
#define ssygv_f77 F77_FUNC(ssygv, SSYGV)
#define ssygvd_f77 F77_FUNC(ssygvd, SSYGVD)
#define ssygvx_f77 F77_FUNC(ssygvx, SSYGVX)
#define ssyrfs_f77 F77_FUNC(ssyrfs, SSYRFS)
#define ssysv_f77 F77_FUNC(ssysv, SSYSV)
#define ssysvx_f77 F77_FUNC(ssysvx, SSYSVX)
#define ssytd2_f77 F77_FUNC(ssytd2, SSYTD2)
#define ssytf2_f77 F77_FUNC(ssytf2, SSYTF2)
#define ssytrd_f77 F77_FUNC(ssytrd, SSYTRD)
#define ssytrf_f77 F77_FUNC(ssytrf, SSYTRF)
#define ssytri_f77 F77_FUNC(ssytri, SSYTRI)
#define ssytrs_f77 F77_FUNC(ssytrs, SSYTRS)
#define stbcon_f77 F77_FUNC(stbcon, STBCON)
#define stbrfs_f77 F77_FUNC(stbrfs, STBRFS)
#define stbtrs_f77 F77_FUNC(stbtrs, STBTRS)
#define stgevc_f77 F77_FUNC(stgevc, STGEVC)
#define stgex2_f77 F77_FUNC(stgex2, STGEX2)
#define stgexc_f77 F77_FUNC(stgexc, STGEXC)
#define stgsen_f77 F77_FUNC(stgsen, STGSEN)
#define stgsja_f77 F77_FUNC(stgsja, STGSJA)
#define stgsna_f77 F77_FUNC(stgsna, STGSNA)
#define stgsy2_f77 F77_FUNC(stgsy2, STGSY2)
#define stgsyl_f77 F77_FUNC(stgsyl, STGSYL)
#define stpcon_f77 F77_FUNC(stpcon, STPCON)
#define stprfs_f77 F77_FUNC(stprfs, STPRFS)
#define stptri_f77 F77_FUNC(stptri, STPTRI)
#define stptrs_f77 F77_FUNC(stptrs, STPTRS)
#define strcon_f77 F77_FUNC(strcon, STRCON)
#define strevc_f77 F77_FUNC(strevc, STREVC)
#define strexc_f77 F77_FUNC(strexc, STREXC)
#define strrfs_f77 F77_FUNC(strrfs, STRRFS)
#define strsen_f77 F77_FUNC(strsen, STRSEN)
#define strsna_f77 F77_FUNC(strsna, STRSNA)
#define strsyl_f77 F77_FUNC(strsyl, STRSYL)
#define strti2_f77 F77_FUNC(strti2, STRTI2)
#define strtri_f77 F77_FUNC(strtri, STRTRI)
#define strtrs_f77 F77_FUNC(strtrs, STRTRS)
#define stzrqf_f77 F77_FUNC(stzrqf, STZRQF)
#define stzrzf_f77 F77_FUNC(stzrzf, STZRZF)
#define xerbla_f77 F77_FUNC(xerbla, XERBLA)
#define zbdsqr_f77 F77_FUNC(zbdsqr, ZBDSQR)
#define zcgesv_f77 F77_FUNC(zcgesv, ZCGESV)
#define zdrscl_f77 F77_FUNC(zdrscl, ZDRSCL)
#define zgbbrd_f77 F77_FUNC(zgbbrd, ZGBBRD)
#define zgbcon_f77 F77_FUNC(zgbcon, ZGBCON)
#define zgbequ_f77 F77_FUNC(zgbequ, ZGBEQU)
#define zgbrfs_f77 F77_FUNC(zgbrfs, ZGBRFS)
#define zgbsv_f77 F77_FUNC(zgbsv, ZGBSV)
#define zgbsvx_f77 F77_FUNC(zgbsvx, ZGBSVX)
#define zgbtf2_f77 F77_FUNC(zgbtf2, ZGBTF2)
#define zgbtrf_f77 F77_FUNC(zgbtrf, ZGBTRF)
#define zgbtrs_f77 F77_FUNC(zgbtrs, ZGBTRS)
#define zgebak_f77 F77_FUNC(zgebak, ZGEBAK)
#define zgebal_f77 F77_FUNC(zgebal, ZGEBAL)
#define zgebd2_f77 F77_FUNC(zgebd2, ZGEBD2)
#define zgebrd_f77 F77_FUNC(zgebrd, ZGEBRD)
#define zgecon_f77 F77_FUNC(zgecon, ZGECON)
#define zgeequ_f77 F77_FUNC(zgeequ, ZGEEQU)
#define zgees_f77 F77_FUNC(zgees, ZGEES)
#define zgeesx_f77 F77_FUNC(zgeesx, ZGEESX)
#define zgeev_f77 F77_FUNC(zgeev, ZGEEV)
#define zgeevx_f77 F77_FUNC(zgeevx, ZGEEVX)
#define zgegs_f77 F77_FUNC(zgegs, ZGEGS)
#define zgegv_f77 F77_FUNC(zgegv, ZGEGV)
#define zgehd2_f77 F77_FUNC(zgehd2, ZGEHD2)
#define zgehrd_f77 F77_FUNC(zgehrd, ZGEHRD)
#define zgelq2_f77 F77_FUNC(zgelq2, ZGELQ2)
#define zgelqf_f77 F77_FUNC(zgelqf, ZGELQF)
#define zgels_f77 F77_FUNC(zgels, ZGELS)
#define zgelsd_f77 F77_FUNC(zgelsd, ZGELSD)
#define zgelss_f77 F77_FUNC(zgelss, ZGELSS)
#define zgelsx_f77 F77_FUNC(zgelsx, ZGELSX)
#define zgelsy_f77 F77_FUNC(zgelsy, ZGELSY)
#define zgeql2_f77 F77_FUNC(zgeql2, ZGEQL2)
#define zgeqlf_f77 F77_FUNC(zgeqlf, ZGEQLF)
#define zgeqp3_f77 F77_FUNC(zgeqp3, ZGEQP3)
#define zgeqpf_f77 F77_FUNC(zgeqpf, ZGEQPF)
#define zgeqr2_f77 F77_FUNC(zgeqr2, ZGEQR2)
#define zgeqrf_f77 F77_FUNC(zgeqrf, ZGEQRF)
#define zgerfs_f77 F77_FUNC(zgerfs, ZGERFS)
#define zgerq2_f77 F77_FUNC(zgerq2, ZGERQ2)
#define zgerqf_f77 F77_FUNC(zgerqf, ZGERQF)
#define zgesc2_f77 F77_FUNC(zgesc2, ZGESC2)
#define zgesdd_f77 F77_FUNC(zgesdd, ZGESDD)
#define zgesv_f77 F77_FUNC(zgesv, ZGESV)
#define zgesvd_f77 F77_FUNC(zgesvd, ZGESVD)
#define zgesvx_f77 F77_FUNC(zgesvx, ZGESVX)
#define zgetc2_f77 F77_FUNC(zgetc2, ZGETC2)
#define zgetf2_f77 F77_FUNC(zgetf2, ZGETF2)
#define zgetrf_f77 F77_FUNC(zgetrf, ZGETRF)
#define zgetri_f77 F77_FUNC(zgetri, ZGETRI)
#define zgetrs_f77 F77_FUNC(zgetrs, ZGETRS)
#define zggbak_f77 F77_FUNC(zggbak, ZGGBAK)
#define zggbal_f77 F77_FUNC(zggbal, ZGGBAL)
#define zgges_f77 F77_FUNC(zgges, ZGGES)
#define zggesx_f77 F77_FUNC(zggesx, ZGGESX)
#define zggev_f77 F77_FUNC(zggev, ZGGEV)
#define zggevx_f77 F77_FUNC(zggevx, ZGGEVX)
#define zggglm_f77 F77_FUNC(zggglm, ZGGGLM)
#define zgghrd_f77 F77_FUNC(zgghrd, ZGGHRD)
#define zgglse_f77 F77_FUNC(zgglse, ZGGLSE)
#define zggqrf_f77 F77_FUNC(zggqrf, ZGGQRF)
#define zggrqf_f77 F77_FUNC(zggrqf, ZGGRQF)
#define zggsvd_f77 F77_FUNC(zggsvd, ZGGSVD)
#define zggsvp_f77 F77_FUNC(zggsvp, ZGGSVP)
#define zgtcon_f77 F77_FUNC(zgtcon, ZGTCON)
#define zgtrfs_f77 F77_FUNC(zgtrfs, ZGTRFS)
#define zgtsv_f77 F77_FUNC(zgtsv, ZGTSV)
#define zgtsvx_f77 F77_FUNC(zgtsvx, ZGTSVX)
#define zgttrf_f77 F77_FUNC(zgttrf, ZGTTRF)
#define zgttrs_f77 F77_FUNC(zgttrs, ZGTTRS)
#define zgtts2_f77 F77_FUNC(zgtts2, ZGTTS2)
#define zhbev_f77 F77_FUNC(zhbev, ZHBEV)
#define zhbevd_f77 F77_FUNC(zhbevd, ZHBEVD)
#define zhbevx_f77 F77_FUNC(zhbevx, ZHBEVX)
#define zhbgst_f77 F77_FUNC(zhbgst, ZHBGST)
#define zhbgv_f77 F77_FUNC(zhbgv, ZHBGV)
#define zhbgvd_f77 F77_FUNC(zhbgvd, ZHBGVD)
#define zhbgvx_f77 F77_FUNC(zhbgvx, ZHBGVX)
#define zhbtrd_f77 F77_FUNC(zhbtrd, ZHBTRD)
#define zhecon_f77 F77_FUNC(zhecon, ZHECON)
#define zheev_f77 F77_FUNC(zheev, ZHEEV)
#define zheevd_f77 F77_FUNC(zheevd, ZHEEVD)
#define zheevr_f77 F77_FUNC(zheevr, ZHEEVR)
#define zheevx_f77 F77_FUNC(zheevx, ZHEEVX)
#define zhegs2_f77 F77_FUNC(zhegs2, ZHEGS2)
#define zhegst_f77 F77_FUNC(zhegst, ZHEGST)
#define zhegv_f77 F77_FUNC(zhegv, ZHEGV)
#define zhegvd_f77 F77_FUNC(zhegvd, ZHEGVD)
#define zhegvx_f77 F77_FUNC(zhegvx, ZHEGVX)
#define zherfs_f77 F77_FUNC(zherfs, ZHERFS)
#define zhesv_f77 F77_FUNC(zhesv, ZHESV)
#define zhesvx_f77 F77_FUNC(zhesvx, ZHESVX)
#define zhetd2_f77 F77_FUNC(zhetd2, ZHETD2)
#define zhetf2_f77 F77_FUNC(zhetf2, ZHETF2)
#define zhetrd_f77 F77_FUNC(zhetrd, ZHETRD)
#define zhetrf_f77 F77_FUNC(zhetrf, ZHETRF)
#define zhetri_f77 F77_FUNC(zhetri, ZHETRI)
#define zhetrs_f77 F77_FUNC(zhetrs, ZHETRS)
#define zhgeqz_f77 F77_FUNC(zhgeqz, ZHGEQZ)
#define zhpcon_f77 F77_FUNC(zhpcon, ZHPCON)
#define zhpev_f77 F77_FUNC(zhpev, ZHPEV)
#define zhpevd_f77 F77_FUNC(zhpevd, ZHPEVD)
#define zhpevx_f77 F77_FUNC(zhpevx, ZHPEVX)
#define zhpgst_f77 F77_FUNC(zhpgst, ZHPGST)
#define zhpgv_f77 F77_FUNC(zhpgv, ZHPGV)
#define zhpgvd_f77 F77_FUNC(zhpgvd, ZHPGVD)
#define zhpgvx_f77 F77_FUNC(zhpgvx, ZHPGVX)
#define zhprfs_f77 F77_FUNC(zhprfs, ZHPRFS)
#define zhpsv_f77 F77_FUNC(zhpsv, ZHPSV)
#define zhpsvx_f77 F77_FUNC(zhpsvx, ZHPSVX)
#define zhptrd_f77 F77_FUNC(zhptrd, ZHPTRD)
#define zhptrf_f77 F77_FUNC(zhptrf, ZHPTRF)
#define zhptri_f77 F77_FUNC(zhptri, ZHPTRI)
#define zhptrs_f77 F77_FUNC(zhptrs, ZHPTRS)
#define zhsein_f77 F77_FUNC(zhsein, ZHSEIN)
#define zhseqr_f77 F77_FUNC(zhseqr, ZHSEQR)
#define zlabrd_f77 F77_FUNC(zlabrd, ZLABRD)
#define zlacgv_f77 F77_FUNC(zlacgv, ZLACGV)
#define zlacn2_f77 F77_FUNC(zlacn2, ZLACN2)
#define zlacon_f77 F77_FUNC(zlacon, ZLACON)
#define zlacp2_f77 F77_FUNC(zlacp2, ZLACP2)
#define zlacpy_f77 F77_FUNC(zlacpy, ZLACPY)
#define zlacrm_f77 F77_FUNC(zlacrm, ZLACRM)
#define zlacrt_f77 F77_FUNC(zlacrt, ZLACRT)
#define zladiv_f77 F77_FUNC(zladiv, ZLADIV)
#define zlaed0_f77 F77_FUNC(zlaed0, ZLAED0)
#define zlaed7_f77 F77_FUNC(zlaed7, ZLAED7)
#define zlaed8_f77 F77_FUNC(zlaed8, ZLAED8)
#define zlaein_f77 F77_FUNC(zlaein, ZLAEIN)
#define zlaesy_f77 F77_FUNC(zlaesy, ZLAESY)
#define zlaev2_f77 F77_FUNC(zlaev2, ZLAEV2)
#define zlag2c_f77 F77_FUNC(zlag2c, ZLAG2C)
#define zlags2_f77 F77_FUNC(zlags2, ZLAGS2)
#define zlagtm_f77 F77_FUNC(zlagtm, ZLAGTM)
#define zlahef_f77 F77_FUNC(zlahef, ZLAHEF)
#define zlahqr_f77 F77_FUNC(zlahqr, ZLAHQR)
#define zlahr2_f77 F77_FUNC(zlahr2, ZLAHR2)
#define zlahrd_f77 F77_FUNC(zlahrd, ZLAHRD)
#define zlaic1_f77 F77_FUNC(zlaic1, ZLAIC1)
#define zlals0_f77 F77_FUNC(zlals0, ZLALS0)
#define zlalsa_f77 F77_FUNC(zlalsa, ZLALSA)
#define zlalsd_f77 F77_FUNC(zlalsd, ZLALSD)
#define zlangb_f77 F77_FUNC(zlangb, ZLANGB)
#define zlange_f77 F77_FUNC(zlange, ZLANGE)
#define zlangt_f77 F77_FUNC(zlangt, ZLANGT)
#define zlanhb_f77 F77_FUNC(zlanhb, ZLANHB)
#define zlanhe_f77 F77_FUNC(zlanhe, ZLANHE)
#define zlanhp_f77 F77_FUNC(zlanhp, ZLANHP)
#define zlanhs_f77 F77_FUNC(zlanhs, ZLANHS)
#define zlanht_f77 F77_FUNC(zlanht, ZLANHT)
#define zlansb_f77 F77_FUNC(zlansb, ZLANSB)
#define zlansp_f77 F77_FUNC(zlansp, ZLANSP)
#define zlansy_f77 F77_FUNC(zlansy, ZLANSY)
#define zlantb_f77 F77_FUNC(zlantb, ZLANTB)
#define zlantp_f77 F77_FUNC(zlantp, ZLANTP)
#define zlantr_f77 F77_FUNC(zlantr, ZLANTR)
#define zlapll_f77 F77_FUNC(zlapll, ZLAPLL)
#define zlapmt_f77 F77_FUNC(zlapmt, ZLAPMT)
#define zlaqgb_f77 F77_FUNC(zlaqgb, ZLAQGB)
#define zlaqge_f77 F77_FUNC(zlaqge, ZLAQGE)
#define zlaqhb_f77 F77_FUNC(zlaqhb, ZLAQHB)
#define zlaqhe_f77 F77_FUNC(zlaqhe, ZLAQHE)
#define zlaqhp_f77 F77_FUNC(zlaqhp, ZLAQHP)
#define zlaqp2_f77 F77_FUNC(zlaqp2, ZLAQP2)
#define zlaqps_f77 F77_FUNC(zlaqps, ZLAQPS)
#define zlaqr0_f77 F77_FUNC(zlaqr0, ZLAQR0)
#define zlaqr1_f77 F77_FUNC(zlaqr1, ZLAQR1)
#define zlaqr2_f77 F77_FUNC(zlaqr2, ZLAQR2)
#define zlaqr3_f77 F77_FUNC(zlaqr3, ZLAQR3)
#define zlaqr4_f77 F77_FUNC(zlaqr4, ZLAQR4)
#define zlaqr5_f77 F77_FUNC(zlaqr5, ZLAQR5)
#define zlaqsb_f77 F77_FUNC(zlaqsb, ZLAQSB)
#define zlaqsp_f77 F77_FUNC(zlaqsp, ZLAQSP)
#define zlaqsy_f77 F77_FUNC(zlaqsy, ZLAQSY)
#define zlar1v_f77 F77_FUNC(zlar1v, ZLAR1V)
#define zlar2v_f77 F77_FUNC(zlar2v, ZLAR2V)
#define zlarcm_f77 F77_FUNC(zlarcm, ZLARCM)
#define zlarf_f77 F77_FUNC(zlarf, ZLARF)
#define zlarfb_f77 F77_FUNC(zlarfb, ZLARFB)
#define zlarfg_f77 F77_FUNC(zlarfg, ZLARFG)
#define zlarft_f77 F77_FUNC(zlarft, ZLARFT)
#define zlarfx_f77 F77_FUNC(zlarfx, ZLARFX)
#define zlargv_f77 F77_FUNC(zlargv, ZLARGV)
#define zlarnv_f77 F77_FUNC(zlarnv, ZLARNV)
#define zlarrv_f77 F77_FUNC(zlarrv, ZLARRV)
#define zlartg_f77 F77_FUNC(zlartg, ZLARTG)
#define zlartv_f77 F77_FUNC(zlartv, ZLARTV)
#define zlarz_f77 F77_FUNC(zlarz, ZLARZ)
#define zlarzb_f77 F77_FUNC(zlarzb, ZLARZB)
#define zlarzt_f77 F77_FUNC(zlarzt, ZLARZT)
#define zlascl_f77 F77_FUNC(zlascl, ZLASCL)
#define zlaset_f77 F77_FUNC(zlaset, ZLASET)
#define zlasr_f77 F77_FUNC(zlasr, ZLASR)
#define zlassq_f77 F77_FUNC(zlassq, ZLASSQ)
#define zlaswp_f77 F77_FUNC(zlaswp, ZLASWP)
#define zlasyf_f77 F77_FUNC(zlasyf, ZLASYF)
#define zlatbs_f77 F77_FUNC(zlatbs, ZLATBS)
#define zlatdf_f77 F77_FUNC(zlatdf, ZLATDF)
#define zlatps_f77 F77_FUNC(zlatps, ZLATPS)
#define zlatrd_f77 F77_FUNC(zlatrd, ZLATRD)
#define zlatrs_f77 F77_FUNC(zlatrs, ZLATRS)
#define zlatrz_f77 F77_FUNC(zlatrz, ZLATRZ)
#define zlatzm_f77 F77_FUNC(zlatzm, ZLATZM)
#define zlauu2_f77 F77_FUNC(zlauu2, ZLAUU2)
#define zlauum_f77 F77_FUNC(zlauum, ZLAUUM)
#define zpbcon_f77 F77_FUNC(zpbcon, ZPBCON)
#define zpbequ_f77 F77_FUNC(zpbequ, ZPBEQU)
#define zpbrfs_f77 F77_FUNC(zpbrfs, ZPBRFS)
#define zpbstf_f77 F77_FUNC(zpbstf, ZPBSTF)
#define zpbsv_f77 F77_FUNC(zpbsv, ZPBSV)
#define zpbsvx_f77 F77_FUNC(zpbsvx, ZPBSVX)
#define zpbtf2_f77 F77_FUNC(zpbtf2, ZPBTF2)
#define zpbtrf_f77 F77_FUNC(zpbtrf, ZPBTRF)
#define zpbtrs_f77 F77_FUNC(zpbtrs, ZPBTRS)
#define zpocon_f77 F77_FUNC(zpocon, ZPOCON)
#define zpoequ_f77 F77_FUNC(zpoequ, ZPOEQU)
#define zporfs_f77 F77_FUNC(zporfs, ZPORFS)
#define zposv_f77 F77_FUNC(zposv, ZPOSV)
#define zposvx_f77 F77_FUNC(zposvx, ZPOSVX)
#define zpotf2_f77 F77_FUNC(zpotf2, ZPOTF2)
#define zpotrf_f77 F77_FUNC(zpotrf, ZPOTRF)
#define zpotri_f77 F77_FUNC(zpotri, ZPOTRI)
#define zpotrs_f77 F77_FUNC(zpotrs, ZPOTRS)
#define zppcon_f77 F77_FUNC(zppcon, ZPPCON)
#define zppequ_f77 F77_FUNC(zppequ, ZPPEQU)
#define zpprfs_f77 F77_FUNC(zpprfs, ZPPRFS)
#define zppsv_f77 F77_FUNC(zppsv, ZPPSV)
#define zppsvx_f77 F77_FUNC(zppsvx, ZPPSVX)
#define zpptrf_f77 F77_FUNC(zpptrf, ZPPTRF)
#define zpptri_f77 F77_FUNC(zpptri, ZPPTRI)
#define zpptrs_f77 F77_FUNC(zpptrs, ZPPTRS)
#define zptcon_f77 F77_FUNC(zptcon, ZPTCON)
#define zpteqr_f77 F77_FUNC(zpteqr, ZPTEQR)
#define zptrfs_f77 F77_FUNC(zptrfs, ZPTRFS)
#define zptsv_f77 F77_FUNC(zptsv, ZPTSV)
#define zptsvx_f77 F77_FUNC(zptsvx, ZPTSVX)
#define zpttrf_f77 F77_FUNC(zpttrf, ZPTTRF)
#define zpttrs_f77 F77_FUNC(zpttrs, ZPTTRS)
#define zptts2_f77 F77_FUNC(zptts2, ZPTTS2)
#define zrot_f77 F77_FUNC(zrot, ZROT)
#define zspcon_f77 F77_FUNC(zspcon, ZSPCON)
#define zspmv_f77 F77_FUNC(zspmv, ZSPMV)
#define zspr_f77 F77_FUNC(zspr, ZSPR)
#define zsprfs_f77 F77_FUNC(zsprfs, ZSPRFS)
#define zspsv_f77 F77_FUNC(zspsv, ZSPSV)
#define zspsvx_f77 F77_FUNC(zspsvx, ZSPSVX)
#define zsptrf_f77 F77_FUNC(zsptrf, ZSPTRF)
#define zsptri_f77 F77_FUNC(zsptri, ZSPTRI)
#define zsptrs_f77 F77_FUNC(zsptrs, ZSPTRS)
#define zstedc_f77 F77_FUNC(zstedc, ZSTEDC)
#define zstegr_f77 F77_FUNC(zstegr, ZSTEGR)
#define zstein_f77 F77_FUNC(zstein, ZSTEIN)
#define zstemr_f77 F77_FUNC(zstemr, ZSTEMR)
#define zsteqr_f77 F77_FUNC(zsteqr, ZSTEQR)
#define zsycon_f77 F77_FUNC(zsycon, ZSYCON)
#define zsymv_f77 F77_FUNC(zsymv, ZSYMV)
#define zsyr_f77 F77_FUNC(zsyr, ZSYR)
#define zsyrfs_f77 F77_FUNC(zsyrfs, ZSYRFS)
#define zsysv_f77 F77_FUNC(zsysv, ZSYSV)
#define zsysvx_f77 F77_FUNC(zsysvx, ZSYSVX)
#define zsytf2_f77 F77_FUNC(zsytf2, ZSYTF2)
#define zsytrf_f77 F77_FUNC(zsytrf, ZSYTRF)
#define zsytri_f77 F77_FUNC(zsytri, ZSYTRI)
#define zsytrs_f77 F77_FUNC(zsytrs, ZSYTRS)
#define ztbcon_f77 F77_FUNC(ztbcon, ZTBCON)
#define ztbrfs_f77 F77_FUNC(ztbrfs, ZTBRFS)
#define ztbtrs_f77 F77_FUNC(ztbtrs, ZTBTRS)
#define ztgevc_f77 F77_FUNC(ztgevc, ZTGEVC)
#define ztgex2_f77 F77_FUNC(ztgex2, ZTGEX2)
#define ztgexc_f77 F77_FUNC(ztgexc, ZTGEXC)
#define ztgsen_f77 F77_FUNC(ztgsen, ZTGSEN)
#define ztgsja_f77 F77_FUNC(ztgsja, ZTGSJA)
#define ztgsna_f77 F77_FUNC(ztgsna, ZTGSNA)
#define ztgsy2_f77 F77_FUNC(ztgsy2, ZTGSY2)
#define ztgsyl_f77 F77_FUNC(ztgsyl, ZTGSYL)
#define ztpcon_f77 F77_FUNC(ztpcon, ZTPCON)
#define ztprfs_f77 F77_FUNC(ztprfs, ZTPRFS)
#define ztptri_f77 F77_FUNC(ztptri, ZTPTRI)
#define ztptrs_f77 F77_FUNC(ztptrs, ZTPTRS)
#define ztrcon_f77 F77_FUNC(ztrcon, ZTRCON)
#define ztrevc_f77 F77_FUNC(ztrevc, ZTREVC)
#define ztrexc_f77 F77_FUNC(ztrexc, ZTREXC)
#define ztrrfs_f77 F77_FUNC(ztrrfs, ZTRRFS)
#define ztrsen_f77 F77_FUNC(ztrsen, ZTRSEN)
#define ztrsna_f77 F77_FUNC(ztrsna, ZTRSNA)
#define ztrsyl_f77 F77_FUNC(ztrsyl, ZTRSYL)
#define ztrti2_f77 F77_FUNC(ztrti2, ZTRTI2)
#define ztrtri_f77 F77_FUNC(ztrtri, ZTRTRI)
#define ztrtrs_f77 F77_FUNC(ztrtrs, ZTRTRS)
#define ztzrqf_f77 F77_FUNC(ztzrqf, ZTZRQF)
#define ztzrzf_f77 F77_FUNC(ztzrzf, ZTZRZF)
#define zung2l_f77 F77_FUNC(zung2l, ZUNG2L)
#define zung2r_f77 F77_FUNC(zung2r, ZUNG2R)
#define zungbr_f77 F77_FUNC(zungbr, ZUNGBR)
#define zunghr_f77 F77_FUNC(zunghr, ZUNGHR)
#define zungl2_f77 F77_FUNC(zungl2, ZUNGL2)
#define zunglq_f77 F77_FUNC(zunglq, ZUNGLQ)
#define zungql_f77 F77_FUNC(zungql, ZUNGQL)
#define zungqr_f77 F77_FUNC(zungqr, ZUNGQR)
#define zungr2_f77 F77_FUNC(zungr2, ZUNGR2)
#define zungrq_f77 F77_FUNC(zungrq, ZUNGRQ)
#define zungtr_f77 F77_FUNC(zungtr, ZUNGTR)
#define zunm2l_f77 F77_FUNC(zunm2l, ZUNM2L)
#define zunm2r_f77 F77_FUNC(zunm2r, ZUNM2R)
#define zunmbr_f77 F77_FUNC(zunmbr, ZUNMBR)
#define zunmhr_f77 F77_FUNC(zunmhr, ZUNMHR)
#define zunml2_f77 F77_FUNC(zunml2, ZUNML2)
#define zunmlq_f77 F77_FUNC(zunmlq, ZUNMLQ)
#define zunmql_f77 F77_FUNC(zunmql, ZUNMQL)
#define zunmqr_f77 F77_FUNC(zunmqr, ZUNMQR)
#define zunmr2_f77 F77_FUNC(zunmr2, ZUNMR2)
#define zunmr3_f77 F77_FUNC(zunmr3, ZUNMR3)
#define zunmrq_f77 F77_FUNC(zunmrq, ZUNMRQ)
#define zunmrz_f77 F77_FUNC(zunmrz, ZUNMRZ)
#define zunmtr_f77 F77_FUNC(zunmtr, ZUNMTR)
#define zupgtr_f77 F77_FUNC(zupgtr, ZUPGTR)
#define zupmtr_f77 F77_FUNC(zupmtr, ZUPMTR)
F77_RET_I cbdsqr_f77(const char *uplo, int *n, int *ncvt, int *nru, int *ncc, float *d, float *e, singlecomplex *vt, int *ldvt, singlecomplex *u, int *ldu, singlecomplex *c, int *ldc, float *rwork, int *info);

F77_RET_I cgbbrd_f77(const char *vect, int *m, int *n, int *ncc, int *kl, int *ku, singlecomplex *ab, int *ldab, float *d, float *e, singlecomplex *q, int *ldq, singlecomplex *pt, int *ldpt, singlecomplex *c, int *ldc, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgbcon_f77(const char *norm, int *n, int *kl, int *ku, singlecomplex *ab, int *ldab, int *ipiv, float *anorm, float *rcond, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgbequ_f77(int *m, int *n, int *kl, int *ku, singlecomplex *ab, int *ldab, float *r, float *c, float *rowcnd, float *colcnd, float *amax, int *info);

F77_RET_I cgbrfs_f77(const char *trans, int *n, int *kl, int *ku, int *nrhs, singlecomplex *ab, int *ldab, singlecomplex *afb, int *ldafb, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgbsv_f77(int *n, int *kl, int *ku, int *nrhs, singlecomplex *ab, int *ldab, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I cgbsvx_f77(const char *fact, const char *trans, int *n, int *kl, int *ku, int *nrhs, singlecomplex *ab, int *ldab, singlecomplex *afb, int *ldafb, int *ipiv, const char *equed, float *r, float *c, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgbtf2_f77(int *m, int *n, int *kl, int *ku, singlecomplex *ab, int *ldab, int *ipiv, int *info);

F77_RET_I cgbtrf_f77(int *m, int *n, int *kl, int *ku, singlecomplex *ab, int *ldab, int *ipiv, int *info);

F77_RET_I cgbtrs_f77(const char *trans, int *n, int *kl, int *ku, int *nrhs, singlecomplex *ab, int *ldab, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I cgebak_f77(const char *job, const char *side, int *n, int *ilo, int *ihi, float *scale, int *m, singlecomplex *v, int *ldv, int *info);

F77_RET_I cgebal_f77(const char *job, int *n, singlecomplex *a, int *lda, int *ilo, int *ihi, float *scale, int *info);

F77_RET_I cgebd2_f77(int *m, int *n, singlecomplex *a, int *lda, float *d, float *e, singlecomplex *tauq, singlecomplex *taup, singlecomplex *work, int *info);

F77_RET_I cgebrd_f77(int *m, int *n, singlecomplex *a, int *lda, float *d, float *e, singlecomplex *tauq, singlecomplex *taup, singlecomplex *work, int *lwork, int *info);

F77_RET_I cgecon_f77(const char *norm, int *n, singlecomplex *a, int *lda, float *anorm, float *rcond, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgeequ_f77(int *m, int *n, singlecomplex *a, int *lda, float *r, float *c, float *rowcnd, float *colcnd, float *amax, int *info);

F77_RET_I cgees_f77(const char *jobvs, const char *sort, L_fp select, int *n, singlecomplex *a, int *lda, int *sdim, singlecomplex *w, singlecomplex *vs, int *ldvs, singlecomplex *work, int *lwork, float *rwork, logical *bwork, int *info);

F77_RET_I cgeesx_f77(const char *jobvs, const char *sort, L_fp select, const char *sense, int *n, singlecomplex *a, int *lda, int *sdim, singlecomplex *w, singlecomplex *vs, int *ldvs, float *rconde, float *rcondv, singlecomplex *work, int *lwork, float *rwork, logical *bwork, int *info);

F77_RET_I cgeev_f77(const char *jobvl, const char *jobvr, int *n, singlecomplex *a, int *lda, singlecomplex *w, singlecomplex *vl, int *ldvl, singlecomplex *vr, int *ldvr, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I cgeevx_f77(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, int *n, singlecomplex *a, int *lda, singlecomplex *w, singlecomplex *vl, int *ldvl, singlecomplex *vr, int *ldvr, int *ilo, int *ihi, float *scale, float *abnrm, float *rconde, float *rcondv, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I cgegs_f77(const char *jobvsl, const char *jobvsr, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *alpha, singlecomplex *beta, singlecomplex *vsl, int *ldvsl, singlecomplex *vsr, int *ldvsr, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I cgegv_f77(const char *jobvl, const char *jobvr, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *alpha, singlecomplex *beta, singlecomplex *vl, int *ldvl, singlecomplex *vr, int *ldvr, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I cgehd2_f77(int *n, int *ilo, int *ihi, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *info);

F77_RET_I cgehrd_f77(int *n, int *ilo, int *ihi, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cgelq2_f77(int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *info);

F77_RET_I cgelqf_f77(int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cgels_f77(const char *trans, int *m, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *work, int *lwork, int *info);

F77_RET_I cgelsd_f77(int *m, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, float *s, float *rcond, int *rank, singlecomplex *work, int *lwork, float *rwork, int *iwork, int *info);

F77_RET_I cgelss_f77(int *m, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, float *s, float *rcond, int *rank, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I cgelsx_f77(int *m, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, int *jpvt, float *rcond, int *rank, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgelsy_f77(int *m, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, int *jpvt, float *rcond, int *rank, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I cgeql2_f77(int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *info);

F77_RET_I cgeqlf_f77(int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cgeqp3_f77(int *m, int *n, singlecomplex *a, int *lda, int *jpvt, singlecomplex *tau, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I cgeqpf_f77(int *m, int *n, singlecomplex *a, int *lda, int *jpvt, singlecomplex *tau, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgeqr2_f77(int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *info);

F77_RET_I cgeqrf_f77(int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cgerfs_f77(const char *trans, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *af, int *ldaf, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgerq2_f77(int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *info);

F77_RET_I cgerqf_f77(int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cgesc2_f77(int *n, singlecomplex *a, int *lda, singlecomplex *rhs, int *ipiv, int *jpiv, float *scale);

F77_RET_I cgesdd_f77(const char *jobz, int *m, int *n, singlecomplex *a, int *lda, float *s, singlecomplex *u, int *ldu, singlecomplex *vt, int *ldvt, singlecomplex *work, int *lwork, float *rwork, int *iwork, int *info);

F77_RET_I cgesv_f77(int *n, int *nrhs, singlecomplex *a, int *lda, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I cgesvd_f77(const char *jobu, const char *jobvt, int *m, int *n, singlecomplex *a, int *lda, float *s, singlecomplex *u, int *ldu, singlecomplex *vt, int *ldvt, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I cgesvx_f77(const char *fact, const char *trans, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *af, int *ldaf, int *ipiv, const char *equed, float *r, float *c, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgetc2_f77(int *n, singlecomplex *a, int *lda, int *ipiv, int *jpiv, int *info);

F77_RET_I cgetf2_f77(int *m, int *n, singlecomplex *a, int *lda, int *ipiv, int *info);

F77_RET_I cgetrf_f77(int *m, int *n, singlecomplex *a, int *lda, int *ipiv, int *info);

F77_RET_I cgetri_f77(int *n, singlecomplex *a, int *lda, int *ipiv, singlecomplex *work, int *lwork, int *info);

F77_RET_I cgetrs_f77(const char *trans, int *n, int *nrhs, singlecomplex *a, int *lda, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I cggbak_f77(const char *job, const char *side, int *n, int *ilo, int *ihi, float *lscale, float *rscale, int *m, singlecomplex *v, int *ldv, int *info);

F77_RET_I cggbal_f77(const char *job, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, int *ilo, int *ihi, float *lscale, float *rscale, float *work, int *info);

F77_RET_I cgges_f77(const char *jobvsl, const char *jobvsr, const char *sort, L_fp selctg, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, int *sdim, singlecomplex *alpha, singlecomplex *beta, singlecomplex *vsl, int *ldvsl, singlecomplex *vsr, int *ldvsr, singlecomplex *work, int *lwork, float *rwork, logical *bwork, int *info);

F77_RET_I cggesx_f77(const char *jobvsl, const char *jobvsr, const char *sort, L_fp selctg, const char *sense, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, int *sdim, singlecomplex *alpha, singlecomplex *beta, singlecomplex *vsl, int *ldvsl, singlecomplex *vsr, int *ldvsr, float *rconde, float *rcondv, singlecomplex *work, int *lwork, float *rwork, int *iwork, int *liwork, logical *bwork, int *info);

F77_RET_I cggev_f77(const char *jobvl, const char *jobvr, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *alpha, singlecomplex *beta, singlecomplex *vl, int *ldvl, singlecomplex *vr, int *ldvr, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I cggevx_f77(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *alpha, singlecomplex *beta, singlecomplex *vl, int *ldvl, singlecomplex *vr, int *ldvr, int *ilo, int *ihi, float *lscale, float *rscale, float *abnrm, float *bbnrm, float *rconde, float *rcondv, singlecomplex *work, int *lwork, float *rwork, int *iwork, logical *bwork, int *info);

F77_RET_I cggglm_f77(int *n, int *m, int *p, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *d, singlecomplex *x, singlecomplex *y, singlecomplex *work, int *lwork, int *info);

F77_RET_I cgghrd_f77(const char *compq, const char *compz, int *n, int *ilo, int *ihi, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *q, int *ldq, singlecomplex *z, int *ldz, int *info);

F77_RET_I cgglse_f77(int *m, int *n, int *p, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *c, singlecomplex *d, singlecomplex *x, singlecomplex *work, int *lwork, int *info);

F77_RET_I cggqrf_f77(int *n, int *m, int *p, singlecomplex *a, int *lda, singlecomplex *taua, singlecomplex *b, int *ldb, singlecomplex *taub, singlecomplex *work, int *lwork, int *info);

F77_RET_I cggrqf_f77(int *m, int *p, int *n, singlecomplex *a, int *lda, singlecomplex *taua, singlecomplex *b, int *ldb, singlecomplex *taub, singlecomplex *work, int *lwork, int *info);

F77_RET_I cggsvd_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *n, int *p, int *k, int *l, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, float *alpha, float *beta, singlecomplex *u, int *ldu, singlecomplex *v, int *ldv, singlecomplex *q, int *ldq, singlecomplex *work, float *rwork, int *iwork, int *info);

F77_RET_I cggsvp_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *p, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, float *tola, float *tolb, int *k, int *l, singlecomplex *u, int *ldu, singlecomplex *v, int *ldv, singlecomplex *q, int *ldq, int *iwork, float *rwork, singlecomplex *tau, singlecomplex *work, int *info);

F77_RET_I cgtcon_f77(const char *norm, int *n, singlecomplex *dl, singlecomplex *d, singlecomplex *du, singlecomplex *du2, int *ipiv, float *anorm, float *rcond, singlecomplex *work, int *info);

F77_RET_I cgtrfs_f77(const char *trans, int *n, int *nrhs, singlecomplex *dl, singlecomplex *d, singlecomplex *du, singlecomplex *dlf, singlecomplex *df, singlecomplex *duf, singlecomplex *du2, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgtsv_f77(int *n, int *nrhs, singlecomplex *dl, singlecomplex *d, singlecomplex *du, singlecomplex *b, int *ldb, int *info);

F77_RET_I cgtsvx_f77(const char *fact, const char *trans, int *n, int *nrhs, singlecomplex *dl, singlecomplex *d, singlecomplex *du, singlecomplex *dlf, singlecomplex *df, singlecomplex *duf, singlecomplex *du2, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cgttrf_f77(int *n, singlecomplex *dl, singlecomplex *d, singlecomplex *du, singlecomplex *du2, int *ipiv, int *info);

F77_RET_I cgttrs_f77(const char *trans, int *n, int *nrhs, singlecomplex *dl, singlecomplex *d, singlecomplex *du, singlecomplex *du2, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I cgtts2_f77(int *itrans, int *n, int *nrhs, singlecomplex *dl, singlecomplex *d, singlecomplex *du, singlecomplex *du2, int *ipiv, singlecomplex *b, int *ldb);

F77_RET_I chbev_f77(const char *jobz, const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, float *w, singlecomplex *z, int *ldz, singlecomplex *work, float *rwork, int *info);

F77_RET_I chbevd_f77(const char *jobz, const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, float *w, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I chbevx_f77(const char *jobz, const char *range, const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, singlecomplex *q, int *ldq, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, singlecomplex *z, int *ldz, singlecomplex *work, float *rwork, int *iwork, int *ifail, int *info);

F77_RET_I chbgst_f77(const char *vect, const char *uplo, int *n, int *ka, int *kb, singlecomplex *ab, int *ldab, singlecomplex *bb, int *ldbb, singlecomplex *x, int *ldx, singlecomplex *work, float *rwork, int *info);

F77_RET_I chbgv_f77(const char *jobz, const char *uplo, int *n, int *ka, int *kb, singlecomplex *ab, int *ldab, singlecomplex *bb, int *ldbb, float *w, singlecomplex *z, int *ldz, singlecomplex *work, float *rwork, int *info);

F77_RET_I chbgvd_f77(const char *jobz, const char *uplo, int *n, int *ka, int *kb, singlecomplex *ab, int *ldab, singlecomplex *bb, int *ldbb, float *w, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I chbgvx_f77(const char *jobz, const char *range, const char *uplo, int *n, int *ka, int *kb, singlecomplex *ab, int *ldab, singlecomplex *bb, int *ldbb, singlecomplex *q, int *ldq, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, singlecomplex *z, int *ldz, singlecomplex *work, float *rwork, int *iwork, int *ifail, int *info);

F77_RET_I chbtrd_f77(const char *vect, const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, float *d, float *e, singlecomplex *q, int *ldq, singlecomplex *work, int *info);

F77_RET_I checon_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *ipiv, float *anorm, float *rcond, singlecomplex *work, int *info);

F77_RET_I cheev_f77(const char *jobz, const char *uplo, int *n, singlecomplex *a, int *lda, float *w, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I cheevd_f77(const char *jobz, const char *uplo, int *n, singlecomplex *a, int *lda, float *w, singlecomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I cheevr_f77(const char *jobz, const char *range, const char *uplo, int *n, singlecomplex *a, int *lda, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, singlecomplex *z, int *ldz, int *isuppz, singlecomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I cheevx_f77(const char *jobz, const char *range, const char *uplo, int *n, singlecomplex *a, int *lda, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, float *rwork, int *iwork, int *ifail, int *info);

F77_RET_I chegs2_f77(int *itype, const char *uplo, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, int *info);

F77_RET_I chegst_f77(int *itype, const char *uplo, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, int *info);

F77_RET_I chegv_f77(int *itype, const char *jobz, const char *uplo, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, float *w, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I chegvd_f77(int *itype, const char *jobz, const char *uplo, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, float *w, singlecomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I chegvx_f77(int *itype, const char *jobz, const char *range, const char *uplo, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, float *rwork, int *iwork, int *ifail, int *info);

F77_RET_I cherfs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *af, int *ldaf, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I chesv_f77(const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *work, int *lwork, int *info);

F77_RET_I chesvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *af, int *ldaf, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I chetd2_f77(const char *uplo, int *n, singlecomplex *a, int *lda, float *d, float *e, singlecomplex *tau, int *info);

F77_RET_I chetf2_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *ipiv, int *info);

F77_RET_I chetrd_f77(const char *uplo, int *n, singlecomplex *a, int *lda, float *d, float *e, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I chetrf_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *ipiv, singlecomplex *work, int *lwork, int *info);

F77_RET_I chetri_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *ipiv, singlecomplex *work, int *info);

F77_RET_I chetrs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I chgeqz_f77(const char *job, const char *compq, const char *compz, int *n, int *ilo, int *ihi, singlecomplex *h, int *ldh, singlecomplex *t, int *ldt, singlecomplex *alpha, singlecomplex *beta, singlecomplex *q, int *ldq, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I chpcon_f77(const char *uplo, int *n, singlecomplex *ap, int *ipiv, float *anorm, float *rcond, singlecomplex *work, int *info);

F77_RET_I chpev_f77(const char *jobz, const char *uplo, int *n, singlecomplex *ap, float *w, singlecomplex *z, int *ldz, singlecomplex *work, float *rwork, int *info);

F77_RET_I chpevd_f77(const char *jobz, const char *uplo, int *n, singlecomplex *ap, float *w, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I chpevx_f77(const char *jobz, const char *range, const char *uplo, int *n, singlecomplex *ap, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, singlecomplex *z, int *ldz, singlecomplex *work, float *rwork, int *iwork, int *ifail, int *info);

F77_RET_I chpgst_f77(int *itype, const char *uplo, int *n, singlecomplex *ap, singlecomplex *bp, int *info);

F77_RET_I chpgv_f77(int *itype, const char *jobz, const char *uplo, int *n, singlecomplex *ap, singlecomplex *bp, float *w, singlecomplex *z, int *ldz, singlecomplex *work, float *rwork, int *info);

F77_RET_I chpgvd_f77(int *itype, const char *jobz, const char *uplo, int *n, singlecomplex *ap, singlecomplex *bp, float *w, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I chpgvx_f77(int *itype, const char *jobz, const char *range, const char *uplo, int *n, singlecomplex *ap, singlecomplex *bp, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, singlecomplex *z, int *ldz, singlecomplex *work, float *rwork, int *iwork, int *ifail, int *info);

F77_RET_I chprfs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *ap, singlecomplex *afp, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I chpsv_f77(const char *uplo, int *n, int *nrhs, singlecomplex *ap, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I chpsvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, singlecomplex *ap, singlecomplex *afp, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I chptrd_f77(const char *uplo, int *n, singlecomplex *ap, float *d, float *e, singlecomplex *tau, int *info);

F77_RET_I chptrf_f77(const char *uplo, int *n, singlecomplex *ap, int *ipiv, int *info);

F77_RET_I chptri_f77(const char *uplo, int *n, singlecomplex *ap, int *ipiv, singlecomplex *work, int *info);

F77_RET_I chptrs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *ap, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I chsein_f77(const char *side, const char *eigsrc, const char *initv, logical *select, int *n, singlecomplex *h, int *ldh, singlecomplex *w, singlecomplex *vl, int *ldvl, singlecomplex *vr, int *ldvr, int *mm, int *m, singlecomplex *work, float *rwork, int *ifaill, int *ifailr, int *info);

F77_RET_I chseqr_f77(const char *job, const char *compz, int *n, int *ilo, int *ihi, singlecomplex *h, int *ldh, singlecomplex *w, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, int *info);

F77_RET_I clabrd_f77(int *m, int *n, int *nb, singlecomplex *a, int *lda, float *d, float *e, singlecomplex *tauq, singlecomplex *taup, singlecomplex *x, int *ldx, singlecomplex *y, int *ldy);

F77_RET_I clacgv_f77(int *n, singlecomplex *x, int *incx);

F77_RET_I clacn2_f77(int *n, singlecomplex *v, singlecomplex *x, float *est, int *kase, int *isave);

F77_RET_I clacon_f77(int *n, singlecomplex *v, singlecomplex *x, float *est, int *kase);

F77_RET_I clacp2_f77(const char *uplo, int *m, int *n, float *a, int *lda, singlecomplex *b, int *ldb);

F77_RET_I clacpy_f77(const char *uplo, int *m, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb);

F77_RET_I clacrm_f77(int *m, int *n, singlecomplex *a, int *lda, float *b, int *ldb, singlecomplex *c, int *ldc, float *rwork);

F77_RET_I clacrt_f77(int *n, singlecomplex *cx, int *incx, singlecomplex *cy, int *incy, singlecomplex *c, singlecomplex *s);

F77_RET_C cladiv_f77(singlecomplex *ret_val, singlecomplex *x, singlecomplex *y);

F77_RET_I claed0_f77(int *qsiz, int *n, float *d, float *e, singlecomplex *q, int *ldq, singlecomplex *qstore, int *ldqs, float *rwork, int *iwork, int *info);

F77_RET_I claed7_f77(int *n, int *cutpnt, int *qsiz, int *tlvls, int *curlvl, int *curpbm, float *d, singlecomplex *q, int *ldq, float *rho, int *indxq, float *qstore, int *qptr, int *prmptr, int *perm, int *givptr, int *givcol, float *givnum, singlecomplex *work, float *rwork, int *iwork, int *info);

F77_RET_I claed8_f77(int *k, int *n, int *qsiz, singlecomplex *q, int *ldq, float *d, float *rho, int *cutpnt, float *z, float *dlamda, singlecomplex *q2, int *ldq2, float *w, int *indxp, int *indx, int *indxq, int *perm, int *givptr, int *givcol, float *givnum, int *info);

F77_RET_I claein_f77(logical *rightv, logical *noinit, int *n, singlecomplex *h, int *ldh, singlecomplex *w, singlecomplex *v, singlecomplex *b, int *ldb, float *rwork, float *eps3, float *smlnum, int *info);

F77_RET_I claesy_f77(singlecomplex *a, singlecomplex *b, singlecomplex *c, singlecomplex *rt1, singlecomplex *rt2, singlecomplex *evscal, singlecomplex *cs1, singlecomplex *sn1);

F77_RET_I claev2_f77(singlecomplex *a, singlecomplex *b, singlecomplex *c, float *rt1, float *rt2, float *cs1, singlecomplex *sn1);

F77_RET_I clag2z_f77(int *m, int *n, singlecomplex *sa, int *ldsa, doublecomplex *a, int *lda, int *info);

F77_RET_I clags2_f77(logical *upper, float *a1, singlecomplex *a2, float *a3, float *b1, singlecomplex *b2, float *b3, float *csu, singlecomplex *snu, float *csv, singlecomplex *snv, float *csq, singlecomplex *snq);

F77_RET_I clagtm_f77(const char *trans, int *n, int *nrhs, float *alpha, singlecomplex *dl, singlecomplex *d, singlecomplex *du, singlecomplex *x, int *ldx, float *beta, singlecomplex *b, int *ldb);

F77_RET_I clahef_f77(const char *uplo, int *n, int *nb, int *kb, singlecomplex *a, int *lda, int *ipiv, singlecomplex *w, int *ldw, int *info);

F77_RET_I clahqr_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, singlecomplex *h, int *ldh, singlecomplex *w, int *iloz, int *ihiz, singlecomplex *z, int *ldz, int *info);

F77_RET_I clahr2_f77(int *n, int *k, int *nb, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *t, int *ldt, singlecomplex *y, int *ldy);

F77_RET_I clahrd_f77(int *n, int *k, int *nb, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *t, int *ldt, singlecomplex *y, int *ldy);

F77_RET_I claic1_f77(int *job, int *j, singlecomplex *x, float *sest, singlecomplex *w, singlecomplex *gamma, float *sestpr, singlecomplex *s, singlecomplex *c);

F77_RET_I clals0_f77(int *icompq, int *nl, int *nr, int *sqre, int *nrhs, singlecomplex *b, int *ldb, singlecomplex *bx, int *ldbx, int *perm, int *givptr, int *givcol, int *ldgcol, float *givnum, int *ldgnum, float *poles, float *difl, float *difr, float *z, int *k, float *c, float *s, float *rwork, int *info);

F77_RET_I clalsa_f77(int *icompq, int *smlsiz, int *n, int *nrhs, singlecomplex *b, int *ldb, singlecomplex *bx, int *ldbx, float *u, int *ldu, float *vt, int *k, float *difl, float *difr, float *z, float *poles, int *givptr, int *givcol, int *ldgcol, int *perm, float *givnum, float *c, float *s, float *rwork, int *iwork, int *info);

F77_RET_I clalsd_f77(const char *uplo, int *smlsiz, int *n, int *nrhs, float *d, float *e, singlecomplex *b, int *ldb, float *rcond, int *rank, singlecomplex *work, float *rwork, int *iwork, int *info);

F77_RET_F clangb_f77(const char *norm, int *n, int *kl, int *ku, singlecomplex *ab, int *ldab, float *work);

F77_RET_F clange_f77(const char *norm, int *m, int *n, singlecomplex *a, int *lda, float *work);

F77_RET_F clangt_f77(const char *norm, int *n, singlecomplex *dl, singlecomplex *d, singlecomplex *du);

F77_RET_F clanhb_f77(const char *norm, const char *uplo, int *n, int *k, singlecomplex *ab, int *ldab, float *work);

F77_RET_F clanhe_f77(const char *norm, const char *uplo, int *n, singlecomplex *a, int *lda, float *work);

F77_RET_F clanhp_f77(const char *norm, const char *uplo, int *n, singlecomplex *ap, float *work);

F77_RET_F clanhs_f77(const char *norm, int *n, singlecomplex *a, int *lda, float *work);

F77_RET_F clanht_f77(const char *norm, int *n, float *d, singlecomplex *e);

F77_RET_F clansb_f77(const char *norm, const char *uplo, int *n, int *k, singlecomplex *ab, int *ldab, float *work);

F77_RET_F clansp_f77(const char *norm, const char *uplo, int *n, singlecomplex *ap, float *work);

F77_RET_F clansy_f77(const char *norm, const char *uplo, int *n, singlecomplex *a, int *lda, float *work);

F77_RET_F clantb_f77(const char *norm, const char *uplo, const char *diag, int *n, int *k, singlecomplex *ab, int *ldab, float *work);

F77_RET_F clantp_f77(const char *norm, const char *uplo, const char *diag, int *n, singlecomplex *ap, float *work);

F77_RET_F clantr_f77(const char *norm, const char *uplo, const char *diag, int *m, int *n, singlecomplex *a, int *lda, float *work);

F77_RET_I clapll_f77(int *n, singlecomplex *x, int *incx, singlecomplex *y, int *incy, float *ssmin);

F77_RET_I clapmt_f77(logical *forwrd, int *m, int *n, singlecomplex *x, int *ldx, int *k);

F77_RET_I claqgb_f77(int *m, int *n, int *kl, int *ku, singlecomplex *ab, int *ldab, float *r, float *c, float *rowcnd, float *colcnd, float *amax, const char *equed);

F77_RET_I claqge_f77(int *m, int *n, singlecomplex *a, int *lda, float *r, float *c, float *rowcnd, float *colcnd, float *amax, const char *equed);

F77_RET_I claqhb_f77(const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, float *s, float *scond, float *amax, const char *equed);

F77_RET_I claqhe_f77(const char *uplo, int *n, singlecomplex *a, int *lda, float *s, float *scond, float *amax, const char *equed);

F77_RET_I claqhp_f77(const char *uplo, int *n, singlecomplex *ap, float *s, float *scond, float *amax, const char *equed);

F77_RET_I claqp2_f77(int *m, int *n, int *offset, singlecomplex *a, int *lda, int *jpvt, singlecomplex *tau, float *vn1, float *vn2, singlecomplex *work);

F77_RET_I claqps_f77(int *m, int *n, int *offset, int *nb, int *kb, singlecomplex *a, int *lda, int *jpvt, singlecomplex *tau, float *vn1, float *vn2, singlecomplex *auxv, singlecomplex *f, int *ldf);

F77_RET_I claqr0_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, singlecomplex *h, int *ldh, singlecomplex *w, int *iloz, int *ihiz, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, int *info);

F77_RET_I claqr1_f77(int *n, singlecomplex *h, int *ldh, singlecomplex *s1, singlecomplex *s2, singlecomplex *v);

F77_RET_I claqr2_f77(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, singlecomplex *h, int *ldh, int *iloz, int *ihiz, singlecomplex *z, int *ldz, int *ns, int *nd, singlecomplex *sh, singlecomplex *v, int *ldv, int *nh, singlecomplex *t, int *ldt, int *nv, singlecomplex *wv, int *ldwv, singlecomplex *work, int *lwork);

F77_RET_I claqr3_f77(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, singlecomplex *h, int *ldh, int *iloz, int *ihiz, singlecomplex *z, int *ldz, int *ns, int *nd, singlecomplex *sh, singlecomplex *v, int *ldv, int *nh, singlecomplex *t, int *ldt, int *nv, singlecomplex *wv, int *ldwv, singlecomplex *work, int *lwork);

F77_RET_I claqr4_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, singlecomplex *h, int *ldh, singlecomplex *w, int *iloz, int *ihiz, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, int *info);

F77_RET_I claqr5_f77(logical *wantt, logical *wantz, int *kacc22, int *n, int *ktop, int *kbot, int *nshfts, singlecomplex *s, singlecomplex *h, int *ldh, int *iloz, int *ihiz, singlecomplex *z, int *ldz, singlecomplex *v, int *ldv, singlecomplex *u, int *ldu, int *nv, singlecomplex *wv, int *ldwv, int *nh, singlecomplex *wh, int *ldwh);

F77_RET_I claqsb_f77(const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, float *s, float *scond, float *amax, const char *equed);

F77_RET_I claqsp_f77(const char *uplo, int *n, singlecomplex *ap, float *s, float *scond, float *amax, const char *equed);

F77_RET_I claqsy_f77(const char *uplo, int *n, singlecomplex *a, int *lda, float *s, float *scond, float *amax, const char *equed);

F77_RET_I clar1v_f77(int *n, int *b1, int *bn, float *lambda, float *d, float *l, float *ld, float *lld, float *pivmin, float *gaptol, singlecomplex *z, logical *wantnc, int *negcnt, float *ztz, float *mingma, int *r, int *isuppz, float *nrminv, float *resid, float *rqcorr, float *work);

F77_RET_I clar2v_f77(int *n, singlecomplex *x, singlecomplex *y, singlecomplex *z, int *incx, float *c, singlecomplex *s, int *incc);

F77_RET_I clarcm_f77(int *m, int *n, float *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *c, int *ldc, float *rwork);

F77_RET_I clarf_f77(const char *side, int *m, int *n, singlecomplex *v, int *incv, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work);

F77_RET_I clarfb_f77(const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n, int *k, singlecomplex *v, int *ldv, singlecomplex *t, int *ldt, singlecomplex *c, int *ldc, singlecomplex *work, int *ldwork);

F77_RET_I clarfg_f77(int *n, singlecomplex *alpha, singlecomplex *x, int *incx, singlecomplex *tau);

F77_RET_I clarft_f77(const char *direct, const char *storev, int *n, int *k, singlecomplex *v, int *ldv, singlecomplex *tau, singlecomplex *t, int *ldt);

F77_RET_I clarfx_f77(const char *side, int *m, int *n, singlecomplex *v, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work);

F77_RET_I clargv_f77(int *n, singlecomplex *x, int *incx, singlecomplex *y, int *incy, float *c, int *incc);

F77_RET_I clarnv_f77(int *idist, int *iseed, int *n, singlecomplex *x);

F77_RET_I clarrv_f77(int *n, float *vl, float *vu, float *d, float *l, float *pivmin, int *isplit, int *m, int *dol, int *dou, float *minrgp, float *rtol1, float *rtol2, float *w, float *werr, float *wgap, int *iblock, int *indexw, float *gers, singlecomplex *z, int *ldz, int *isuppz, float *work, int *iwork, int *info);

F77_RET_I clartg_f77(singlecomplex *f, singlecomplex *g, float *cs, singlecomplex *sn, singlecomplex *r);

F77_RET_I clartv_f77(int *n, singlecomplex *x, int *incx, singlecomplex *y, int *incy, float *c, singlecomplex *s, int *incc);

F77_RET_I clarz_f77(const char *side, int *m, int *n, int *l, singlecomplex *v, int *incv, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work);

F77_RET_I clarzb_f77(const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n, int *k, int *l, singlecomplex *v, int *ldv, singlecomplex *t, int *ldt, singlecomplex *c, int *ldc, singlecomplex *work, int *ldwork);

F77_RET_I clarzt_f77(const char *direct, const char *storev, int *n, int *k, singlecomplex *v, int *ldv, singlecomplex *tau, singlecomplex *t, int *ldt);

F77_RET_I clascl_f77(const char *type, int *kl, int *ku, float *cfrom, float *cto, int *m, int *n, singlecomplex *a, int *lda, int *info);

F77_RET_I claset_f77(const char *uplo, int *m, int *n, singlecomplex *alpha, singlecomplex *beta, singlecomplex *a, int *lda);

F77_RET_I clasr_f77(const char *side, const char *pivot, const char *direct, int *m, int *n, float *c, float *s, singlecomplex *a, int *lda);

F77_RET_I classq_f77(int *n, singlecomplex *x, int *incx, float *scale, float *sumsq);

F77_RET_I claswp_f77(int *n, singlecomplex *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);

F77_RET_I clasyf_f77(const char *uplo, int *n, int *nb, int *kb, singlecomplex *a, int *lda, int *ipiv, singlecomplex *w, int *ldw, int *info);

F77_RET_I clatbs_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, int *kd, singlecomplex *ab, int *ldab, singlecomplex *x, float *scale, float *cnorm, int *info);

F77_RET_I clatdf_f77(int *ijob, int *n, singlecomplex *z, int *ldz, singlecomplex *rhs, float *rdsum, float *rdscal, int *ipiv, int *jpiv);

F77_RET_I clatps_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, singlecomplex *ap, singlecomplex *x, float *scale, float *cnorm, int *info);

F77_RET_I clatrd_f77(const char *uplo, int *n, int *nb, singlecomplex *a, int *lda, float *e, singlecomplex *tau, singlecomplex *w, int *ldw);

F77_RET_I clatrs_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, singlecomplex *a, int *lda, singlecomplex *x, float *scale, float *cnorm, int *info);

F77_RET_I clatrz_f77(int *m, int *n, int *l, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work);

F77_RET_I clatzm_f77(const char *side, int *m, int *n, singlecomplex *v, int *incv, singlecomplex *tau, singlecomplex *c1, singlecomplex *c2, int *ldc, singlecomplex *work);

F77_RET_I clauu2_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *info);

F77_RET_I clauum_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *info);

F77_RET_I cpbcon_f77(const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, float *anorm, float *rcond, singlecomplex *work, float *rwork, int *info);

F77_RET_I cpbequ_f77(const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, float *s, float *scond, float *amax, int *info);

F77_RET_I cpbrfs_f77(const char *uplo, int *n, int *kd, int *nrhs, singlecomplex *ab, int *ldab, singlecomplex *afb, int *ldafb, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cpbstf_f77(const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, int *info);

F77_RET_I cpbsv_f77(const char *uplo, int *n, int *kd, int *nrhs, singlecomplex *ab, int *ldab, singlecomplex *b, int *ldb, int *info);

F77_RET_I cpbsvx_f77(const char *fact, const char *uplo, int *n, int *kd, int *nrhs, singlecomplex *ab, int *ldab, singlecomplex *afb, int *ldafb, const char *equed, float *s, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cpbtf2_f77(const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, int *info);

F77_RET_I cpbtrf_f77(const char *uplo, int *n, int *kd, singlecomplex *ab, int *ldab, int *info);

F77_RET_I cpbtrs_f77(const char *uplo, int *n, int *kd, int *nrhs, singlecomplex *ab, int *ldab, singlecomplex *b, int *ldb, int *info);

F77_RET_I cpocon_f77(const char *uplo, int *n, singlecomplex *a, int *lda, float *anorm, float *rcond, singlecomplex *work, float *rwork, int *info);

F77_RET_I cpoequ_f77(int *n, singlecomplex *a, int *lda, float *s, float *scond, float *amax, int *info);

F77_RET_I cporfs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *af, int *ldaf, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cposv_f77(const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, int *info);

F77_RET_I cposvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *af, int *ldaf, const char *equed, float *s, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cpotf2_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *info);

F77_RET_I cpotrf_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *info);

F77_RET_I cpotri_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *info);

F77_RET_I cpotrs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, int *info);

F77_RET_I cppcon_f77(const char *uplo, int *n, singlecomplex *ap, float *anorm, float *rcond, singlecomplex *work, float *rwork, int *info);

F77_RET_I cppequ_f77(const char *uplo, int *n, singlecomplex *ap, float *s, float *scond, float *amax, int *info);

F77_RET_I cpprfs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *ap, singlecomplex *afp, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cppsv_f77(const char *uplo, int *n, int *nrhs, singlecomplex *ap, singlecomplex *b, int *ldb, int *info);

F77_RET_I cppsvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, singlecomplex *ap, singlecomplex *afp, const char *equed, float *s, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cpptrf_f77(const char *uplo, int *n, singlecomplex *ap, int *info);

F77_RET_I cpptri_f77(const char *uplo, int *n, singlecomplex *ap, int *info);

F77_RET_I cpptrs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *ap, singlecomplex *b, int *ldb, int *info);

F77_RET_I cptcon_f77(int *n, float *d, singlecomplex *e, float *anorm, float *rcond, float *rwork, int *info);

F77_RET_I cpteqr_f77(const char *compz, int *n, float *d, float *e, singlecomplex *z, int *ldz, float *work, int *info);

F77_RET_I cptrfs_f77(const char *uplo, int *n, int *nrhs, float *d, singlecomplex *e, float *df, singlecomplex *ef, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cptsv_f77(int *n, int *nrhs, float *d, singlecomplex *e, singlecomplex *b, int *ldb, int *info);

F77_RET_I cptsvx_f77(const char *fact, int *n, int *nrhs, float *d, singlecomplex *e, float *df, singlecomplex *ef, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cpttrf_f77(int *n, float *d, singlecomplex *e, int *info);

F77_RET_I cpttrs_f77(const char *uplo, int *n, int *nrhs, float *d, singlecomplex *e, singlecomplex *b, int *ldb, int *info);

F77_RET_I cptts2_f77(int *iuplo, int *n, int *nrhs, float *d, singlecomplex *e, singlecomplex *b, int *ldb);

F77_RET_I crot_f77(int *n, singlecomplex *cx, int *incx, singlecomplex *cy, int *incy, float *c, singlecomplex *s);

F77_RET_I cspcon_f77(const char *uplo, int *n, singlecomplex *ap, int *ipiv, float *anorm, float *rcond, singlecomplex *work, int *info);

F77_RET_I cspmv_f77(const char *uplo, int *n, singlecomplex *alpha, singlecomplex *ap, singlecomplex *x, int *incx, singlecomplex *beta, singlecomplex *y, int *incy);

F77_RET_I cspr_f77(const char *uplo, int *n, singlecomplex *alpha, singlecomplex *x, int *incx, singlecomplex *ap);

F77_RET_I csprfs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *ap, singlecomplex *afp, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I cspsv_f77(const char *uplo, int *n, int *nrhs, singlecomplex *ap, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I cspsvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, singlecomplex *ap, singlecomplex *afp, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I csptrf_f77(const char *uplo, int *n, singlecomplex *ap, int *ipiv, int *info);

F77_RET_I csptri_f77(const char *uplo, int *n, singlecomplex *ap, int *ipiv, singlecomplex *work, int *info);

F77_RET_I csptrs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *ap, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I csrscl_f77(int *n, float *sa, singlecomplex *sx, int *incx);

F77_RET_I cstedc_f77(const char *compz, int *n, float *d, float *e, singlecomplex *z, int *ldz, singlecomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I cstegr_f77(const char *jobz, const char *range, int *n, float *d, float *e, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, singlecomplex *z, int *ldz, int *isuppz, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I cstein_f77(int *n, float *d, float *e, int *m, float *w, int *iblock, int *isplit, singlecomplex *z, int *ldz, float *work, int *iwork, int *ifail, int *info);

F77_RET_I cstemr_f77(const char *jobz, const char *range, int *n, float *d, float *e, float *vl, float *vu, int *il, int *iu, int *m, float *w, singlecomplex *z, int *ldz, int *nzc, int *isuppz, logical *tryrac, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I csteqr_f77(const char *compz, int *n, float *d, float *e, singlecomplex *z, int *ldz, float *work, int *info);

F77_RET_I csycon_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *ipiv, float *anorm, float *rcond, singlecomplex *work, int *info);

F77_RET_I csymv_f77(const char *uplo, int *n, singlecomplex *alpha, singlecomplex *a, int *lda, singlecomplex *x, int *incx, singlecomplex *beta, singlecomplex *y, int *incy);

F77_RET_I csyr_f77(const char *uplo, int *n, singlecomplex *alpha, singlecomplex *x, int *incx, singlecomplex *a, int *lda);

F77_RET_I csyrfs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *af, int *ldaf, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I csysv_f77(const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *work, int *lwork, int *info);

F77_RET_I csysvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *af, int *ldaf, int *ipiv, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *rcond, float *ferr, float *berr, singlecomplex *work, int *lwork, float *rwork, int *info);

F77_RET_I csytf2_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *ipiv, int *info);

F77_RET_I csytrf_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *ipiv, singlecomplex *work, int *lwork, int *info);

F77_RET_I csytri_f77(const char *uplo, int *n, singlecomplex *a, int *lda, int *ipiv, singlecomplex *work, int *info);

F77_RET_I csytrs_f77(const char *uplo, int *n, int *nrhs, singlecomplex *a, int *lda, int *ipiv, singlecomplex *b, int *ldb, int *info);

F77_RET_I ctbcon_f77(const char *norm, const char *uplo, const char *diag, int *n, int *kd, singlecomplex *ab, int *ldab, float *rcond, singlecomplex *work, float *rwork, int *info);

F77_RET_I ctbrfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *kd, int *nrhs, singlecomplex *ab, int *ldab, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I ctbtrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *kd, int *nrhs, singlecomplex *ab, int *ldab, singlecomplex *b, int *ldb, int *info);

F77_RET_I ctgevc_f77(const char *side, const char *howmny, logical *select, int *n, singlecomplex *s, int *lds, singlecomplex *p, int *ldp, singlecomplex *vl, int *ldvl, singlecomplex *vr, int *ldvr, int *mm, int *m, singlecomplex *work, float *rwork, int *info);

F77_RET_I ctgex2_f77(logical *wantq, logical *wantz, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *q, int *ldq, singlecomplex *z, int *ldz, int *j1, int *info);

F77_RET_I ctgexc_f77(logical *wantq, logical *wantz, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *q, int *ldq, singlecomplex *z, int *ldz, int *ifst, int *ilst, int *info);

F77_RET_I ctgsen_f77(int *ijob, logical *wantq, logical *wantz, logical *select, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *alpha, singlecomplex *beta, singlecomplex *q, int *ldq, singlecomplex *z, int *ldz, int *m, float *pl, float *pr, float *dif, singlecomplex *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I ctgsja_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *p, int *n, int *k, int *l, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, float *tola, float *tolb, float *alpha, float *beta, singlecomplex *u, int *ldu, singlecomplex *v, int *ldv, singlecomplex *q, int *ldq, singlecomplex *work, int *ncycle, int *info);

F77_RET_I ctgsna_f77(const char *job, const char *howmny, logical *select, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *vl, int *ldvl, singlecomplex *vr, int *ldvr, float *s, float *dif, int *mm, int *m, singlecomplex *work, int *lwork, int *iwork, int *info);

F77_RET_I ctgsy2_f77(const char *trans, int *ijob, int *m, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *c, int *ldc, singlecomplex *d, int *ldd, singlecomplex *e, int *lde, singlecomplex *f, int *ldf, float *scale, float *rdsum, float *rdscal, int *info);

F77_RET_I ctgsyl_f77(const char *trans, int *ijob, int *m, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *c, int *ldc, singlecomplex *d, int *ldd, singlecomplex *e, int *lde, singlecomplex *f, int *ldf, float *scale, float *dif, singlecomplex *work, int *lwork, int *iwork, int *info);

F77_RET_I ctpcon_f77(const char *norm, const char *uplo, const char *diag, int *n, singlecomplex *ap, float *rcond, singlecomplex *work, float *rwork, int *info);

F77_RET_I ctprfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, singlecomplex *ap, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I ctptri_f77(const char *uplo, const char *diag, int *n, singlecomplex *ap, int *info);

F77_RET_I ctptrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, singlecomplex *ap, singlecomplex *b, int *ldb, int *info);

F77_RET_I ctrcon_f77(const char *norm, const char *uplo, const char *diag, int *n, singlecomplex *a, int *lda, float *rcond, singlecomplex *work, float *rwork, int *info);

F77_RET_I ctrevc_f77(const char *side, const char *howmny, logical *select, int *n, singlecomplex *t, int *ldt, singlecomplex *vl, int *ldvl, singlecomplex *vr, int *ldvr, int *mm, int *m, singlecomplex *work, float *rwork, int *info);

F77_RET_I ctrexc_f77(const char *compq, int *n, singlecomplex *t, int *ldt, singlecomplex *q, int *ldq, int *ifst, int *ilst, int *info);

F77_RET_I ctrrfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *x, int *ldx, float *ferr, float *berr, singlecomplex *work, float *rwork, int *info);

F77_RET_I ctrsen_f77(const char *job, const char *compq, logical *select, int *n, singlecomplex *t, int *ldt, singlecomplex *q, int *ldq, singlecomplex *w, int *m, float *s, float *sep, singlecomplex *work, int *lwork, int *info);

F77_RET_I ctrsna_f77(const char *job, const char *howmny, logical *select, int *n, singlecomplex *t, int *ldt, singlecomplex *vl, int *ldvl, singlecomplex *vr, int *ldvr, float *s, float *sep, int *mm, int *m, singlecomplex *work, int *ldwork, float *rwork, int *info);

F77_RET_I ctrsyl_f77(const char *trana, const char *tranb, int *isgn, int *m, int *n, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, singlecomplex *c, int *ldc, float *scale, int *info);

F77_RET_I ctrti2_f77(const char *uplo, const char *diag, int *n, singlecomplex *a, int *lda, int *info);

F77_RET_I ctrtri_f77(const char *uplo, const char *diag, int *n, singlecomplex *a, int *lda, int *info);

F77_RET_I ctrtrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, singlecomplex *a, int *lda, singlecomplex *b, int *ldb, int *info);

F77_RET_I ctzrqf_f77(int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, int *info);

F77_RET_I ctzrzf_f77(int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cung2l_f77(int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *info);

F77_RET_I cung2r_f77(int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *info);

F77_RET_I cungbr_f77(const char *vect, int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cunghr_f77(int *n, int *ilo, int *ihi, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cungl2_f77(int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *info);

F77_RET_I cunglq_f77(int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cungql_f77(int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cungqr_f77(int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cungr2_f77(int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *info);

F77_RET_I cungrq_f77(int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cungtr_f77(const char *uplo, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *work, int *lwork, int *info);

F77_RET_I cunm2l_f77(const char *side, const char *trans, int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *info);

F77_RET_I cunm2r_f77(const char *side, const char *trans, int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *info);

F77_RET_I cunmbr_f77(const char *vect, const char *side, const char *trans, int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *lwork, int *info);

F77_RET_I cunmhr_f77(const char *side, const char *trans, int *m, int *n, int *ilo, int *ihi, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *lwork, int *info);

F77_RET_I cunml2_f77(const char *side, const char *trans, int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *info);

F77_RET_I cunmlq_f77(const char *side, const char *trans, int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *lwork, int *info);

F77_RET_I cunmql_f77(const char *side, const char *trans, int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *lwork, int *info);

F77_RET_I cunmqr_f77(const char *side, const char *trans, int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *lwork, int *info);

F77_RET_I cunmr2_f77(const char *side, const char *trans, int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *info);

F77_RET_I cunmr3_f77(const char *side, const char *trans, int *m, int *n, int *k, int *l, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *info);

F77_RET_I cunmrq_f77(const char *side, const char *trans, int *m, int *n, int *k, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *lwork, int *info);

F77_RET_I cunmrz_f77(const char *side, const char *trans, int *m, int *n, int *k, int *l, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *lwork, int *info);

F77_RET_I cunmtr_f77(const char *side, const char *uplo, const char *trans, int *m, int *n, singlecomplex *a, int *lda, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *lwork, int *info);

F77_RET_I cupgtr_f77(const char *uplo, int *n, singlecomplex *ap, singlecomplex *tau, singlecomplex *q, int *ldq, singlecomplex *work, int *info);

F77_RET_I cupmtr_f77(const char *side, const char *uplo, const char *trans, int *m, int *n, singlecomplex *ap, singlecomplex *tau, singlecomplex *c, int *ldc, singlecomplex *work, int *info);

F77_RET_I dbdsdc_f77(const char *uplo, const char *compq, int *n, double *d, double *e, double *u, int *ldu, double *vt, int *ldvt, double *q, int *iq, double *work, int *iwork, int *info);

F77_RET_I dbdsqr_f77(const char *uplo, int *n, int *ncvt, int *nru, int *ncc, double *d, double *e, double *vt, int *ldvt, double *u, int *ldu, double *c, int *ldc, double *work, int *info);

F77_RET_I ddisna_f77(const char *job, int *m, int *n, double *d, double *sep, int *info);

F77_RET_I dgbbrd_f77(const char *vect, int *m, int *n, int *ncc, int *kl, int *ku, double *ab, int *ldab, double *d, double *e, double *q, int *ldq, double *pt, int *ldpt, double *c, int *ldc, double *work, int *info);

F77_RET_I dgbcon_f77(const char *norm, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dgbequ_f77(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info);

F77_RET_I dgbrfs_f77(const char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, int *ipiv, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dgbsv_f77(int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);

F77_RET_I dgbsvx_f77(const char *fact, const char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, int *ipiv, const char *equed, double *r, double *c, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dgbtf2_f77(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info);

F77_RET_I dgbtrf_f77(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info);

F77_RET_I dgbtrs_f77(const char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);

F77_RET_I dgebak_f77(const char *job, const char *side, int *n, int *ilo, int *ihi, double *scale, int *m, double *v, int *ldv, int *info);

F77_RET_I dgebal_f77(const char *job, int *n, double *a, int *lda, int *ilo, int *ihi, double *scale, int *info);

F77_RET_I dgebd2_f77(int *m, int *n, double *a, int *lda, double *d, double *e, double *tauq, double *taup, double *work, int *info);

F77_RET_I dgebrd_f77(int *m, int *n, double *a, int *lda, double *d, double *e, double *tauq, double *taup, double *work, int *lwork, int *info);

F77_RET_I dgecon_f77(const char *norm, int *n, double *a, int *lda, double *anorm, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dgeequ_f77(int *m, int *n, double *a, int *lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info);

F77_RET_I dgees_f77(const char *jobvs, const char *sort, L_fp select, int *n, double *a, int *lda, int *sdim, double *wr, double *wi, double *vs, int *ldvs, double *work, int *lwork, logical *bwork, int *info);

F77_RET_I dgeesx_f77(const char *jobvs, const char *sort, L_fp select, const char *sense, int *n, double *a, int *lda, int *sdim, double *wr, double *wi, double *vs, int *ldvs, double *rconde, double *rcondv, double *work, int *lwork, int *iwork, int *liwork, logical *bwork, int *info);

F77_RET_I dgeev_f77(const char *jobvl, const char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);

F77_RET_I dgeevx_f77(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, int *ilo, int *ihi, double *scale, double *abnrm, double *rconde, double *rcondv, double *work, int *lwork, int *iwork, int *info);

F77_RET_I dgegs_f77(const char *jobvsl, const char *jobvsr, int *n, double *a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double *beta, double *vsl, int *ldvsl, double *vsr, int *ldvsr, double *work, int *lwork, int *info);

F77_RET_I dgegv_f77(const char *jobvl, const char *jobvr, int *n, double *a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double *beta, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);

F77_RET_I dgehd2_f77(int *n, int *ilo, int *ihi, double *a, int *lda, double *tau, double *work, int *info);

F77_RET_I dgehrd_f77(int *n, int *ilo, int *ihi, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dgelq2_f77(int *m, int *n, double *a, int *lda, double *tau, double *work, int *info);

F77_RET_I dgelqf_f77(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dgels_f77(const char *trans, int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *work, int *lwork, int *info);

F77_RET_I dgelsd_f77(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *s, double *rcond, int *rank, double *work, int *lwork, int *iwork, int *info);

F77_RET_I dgelss_f77(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *s, double *rcond, int *rank, double *work, int *lwork, int *info);

F77_RET_I dgelsx_f77(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *jpvt, double *rcond, int *rank, double *work, int *info);

F77_RET_I dgelsy_f77(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *jpvt, double *rcond, int *rank, double *work, int *lwork, int *info);

F77_RET_I dgeql2_f77(int *m, int *n, double *a, int *lda, double *tau, double *work, int *info);

F77_RET_I dgeqlf_f77(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dgeqp3_f77(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);

F77_RET_I dgeqpf_f77(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *info);

F77_RET_I dgeqr2_f77(int *m, int *n, double *a, int *lda, double *tau, double *work, int *info);

F77_RET_I dgeqrf_f77(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dgerfs_f77(const char *trans, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dgerq2_f77(int *m, int *n, double *a, int *lda, double *tau, double *work, int *info);

F77_RET_I dgerqf_f77(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dgesc2_f77(int *n, double *a, int *lda, double *rhs, int *ipiv, int *jpiv, double *scale);

F77_RET_I dgesdd_f77(const char *jobz, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *iwork, int *info);

F77_RET_I dgesv_f77(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

F77_RET_I dgesvd_f77(const char *jobu, const char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);

F77_RET_I dgesvx_f77(const char *fact, const char *trans, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, const char *equed, double *r, double *c, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dgetc2_f77(int *n, double *a, int *lda, int *ipiv, int *jpiv, int *info);

F77_RET_I dgetf2_f77(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

F77_RET_I dgetrf_f77(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

F77_RET_I dgetri_f77(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

F77_RET_I dgetrs_f77(const char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

F77_RET_I dggbak_f77(const char *job, const char *side, int *n, int *ilo, int *ihi, double *lscale, double *rscale, int *m, double *v, int *ldv, int *info);

F77_RET_I dggbal_f77(const char *job, int *n, double *a, int *lda, double *b, int *ldb, int *ilo, int *ihi, double *lscale, double *rscale, double *work, int *info);

F77_RET_I dgges_f77(const char *jobvsl, const char *jobvsr, const char *sort, L_fp selctg, int *n, double *a, int *lda, double *b, int *ldb, int *sdim, double *alphar, double *alphai, double *beta, double *vsl, int *ldvsl, double *vsr, int *ldvsr, double *work, int *lwork, logical *bwork, int *info);

F77_RET_I dggesx_f77(const char *jobvsl, const char *jobvsr, const char *sort, L_fp selctg, const char *sense, int *n, double *a, int *lda, double *b, int *ldb, int *sdim, double *alphar, double *alphai, double *beta, double *vsl, int *ldvsl, double *vsr, int *ldvsr, double *rconde, double *rcondv, double *work, int *lwork, int *iwork, int *liwork, logical *bwork, int *info);

F77_RET_I dggev_f77(const char *jobvl, const char *jobvr, int *n, double *a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double *beta, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);

F77_RET_I dggevx_f77(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, int *n, double *a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double *beta, double *vl, int *ldvl, double *vr, int *ldvr, int *ilo, int *ihi, double *lscale, double *rscale, double *abnrm, double *bbnrm, double *rconde, double *rcondv, double *work, int *lwork, int *iwork, logical *bwork, int *info);

F77_RET_I dggglm_f77(int *n, int *m, int *p, double *a, int *lda, double *b, int *ldb, double *d, double *x, double *y, double *work, int *lwork, int *info);

F77_RET_I dgghrd_f77(const char *compq, const char *compz, int *n, int *ilo, int *ihi, double *a, int *lda, double *b, int *ldb, double *q, int *ldq, double *z, int *ldz, int *info);

F77_RET_I dgglse_f77(int *m, int *n, int *p, double *a, int *lda, double *b, int *ldb, double *c, double *d, double *x, double *work, int *lwork, int *info);

F77_RET_I dggqrf_f77(int *n, int *m, int *p, double *a, int *lda, double *taua, double *b, int *ldb, double *taub, double *work, int *lwork, int *info);

F77_RET_I dggrqf_f77(int *m, int *p, int *n, double *a, int *lda, double *taua, double *b, int *ldb, double *taub, double *work, int *lwork, int *info);

F77_RET_I dggsvd_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *n, int *p, int *k, int *l, double *a, int *lda, double *b, int *ldb, double *alpha, double *beta, double *u, int *ldu, double *v, int *ldv, double *q, int *ldq, double *work, int *iwork, int *info);

F77_RET_I dggsvp_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *p, int *n, double *a, int *lda, double *b, int *ldb, double *tola, double *tolb, int *k, int *l, double *u, int *ldu, double *v, int *ldv, double *q, int *ldq, int *iwork, double *tau, double *work, int *info);

F77_RET_I dgtcon_f77(const char *norm, int *n, double *dl, double *d, double *du, double *du2, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dgtrfs_f77(const char *trans, int *n, int *nrhs, double *dl, double *d, double *du, double *dlf, double *df, double *duf, double *du2, int *ipiv, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dgtsv_f77(int *n, int *nrhs, double *dl, double *d, double *du, double *b, int *ldb, int *info);

F77_RET_I dgtsvx_f77(const char *fact, const char *trans, int *n, int *nrhs, double *dl, double *d, double *du, double *dlf, double *df, double *duf, double *du2, int *ipiv, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dgttrf_f77(int *n, double *dl, double *d, double *du, double *du2, int *ipiv, int *info);

F77_RET_I dgttrs_f77(const char *trans, int *n, int *nrhs, double *dl, double *d, double *du, double *du2, int *ipiv, double *b, int *ldb, int *info);

F77_RET_I dgtts2_f77(int *itrans, int *n, int *nrhs, double *dl, double *d, double *du, double *du2, int *ipiv, double *b, int *ldb);

F77_RET_I dhgeqz_f77(const char *job, const char *compq, const char *compz, int *n, int *ilo, int *ihi, double *h, int *ldh, double *t, int *ldt, double *alphar, double *alphai, double *beta, double *q, int *ldq, double *z, int *ldz, double *work, int *lwork, int *info);

F77_RET_I dhsein_f77(const char *side, const char *eigsrc, const char *initv, logical *select, int *n, double *h, int *ldh, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, int *mm, int *m, double *work, int *ifaill, int *ifailr, int *info);

F77_RET_I dhseqr_f77(const char *job, const char *compz, int *n, int *ilo, int *ihi, double *h, int *ldh, double *wr, double *wi, double *z, int *ldz, double *work, int *lwork, int *info);

logical disnan_f77(double *din);

F77_RET_I dlabad_f77(double *small, double *large);

F77_RET_I dlabrd_f77(int *m, int *n, int *nb, double *a, int *lda, double *d, double *e, double *tauq, double *taup, double *x, int *ldx, double *y, int *ldy);

F77_RET_I dlacn2_f77(int *n, double *v, double *x, int *isgn, double *est, int *kase, int *isave);

F77_RET_I dlacon_f77(int *n, double *v, double *x, int *isgn, double *est, int *kase);

F77_RET_I dlacpy_f77(const char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb);

F77_RET_I dladiv_f77(double *a, double *b, double *c, double *d, double *p, double *q);

F77_RET_I dlae2_f77(double *a, double *b, double *c, double *rt1, double *rt2);

F77_RET_I dlaebz_f77(int *ijob, int *nitmax, int *n, int *mmax, int *minp, int *nbmin, double *abstol, double *reltol, double *pivmin, double *d, double *e, double *e2, int *nval, double *ab, double *c, int *mout, int *nab, double *work, int *iwork, int *info);

F77_RET_I dlaed0_f77(int *icompq, int *qsiz, int *n, double *d, double *e, double *q, int *ldq, double *qstore, int *ldqs, double *work, int *iwork, int *info);

F77_RET_I dlaed1_f77(int *n, double *d, double *q, int *ldq, int *indxq, double *rho, int *cutpnt, double *work, int *iwork, int *info);

F77_RET_I dlaed2_f77(int *k, int *n, int *n1, double *d, double *q, int *ldq, int *indxq, double *rho, double *z, double *dlamda, double *w, double *q2, int *indx, int *indxc, int *indxp, int *coltyp, int *info);

F77_RET_I dlaed3_f77(int *k, int *n, int *n1, double *d, double *q, int *ldq, double *rho, double *dlamda, double *q2, int *indx, int *ctot, double *w, double *s, int *info);

F77_RET_I dlaed4_f77(int *n, int *i, double *d, double *z, double *delta, double *rho, double *dlam, int *info);

F77_RET_I dlaed5_f77(int *i, double *d, double *z, double *delta, double *rho, double *dlam);

F77_RET_I dlaed6_f77(int *kniter, logical *orgati, double *rho, double *d, double *z, double *finit, double *tau, int *info);

F77_RET_I dlaed7_f77(int *icompq, int *n, int *qsiz, int *tlvls, int *curlvl, int *curpbm, double *d, double *q, int *ldq, int *indxq, double *rho, int *cutpnt, double *qstore, int *qptr, int *prmptr, int *perm, int *givptr, int *givcol, double *givnum, double *work, int *iwork, int *info);

F77_RET_I dlaed8_f77(int *icompq, int *k, int *n, int *qsiz, double *d, double *q, int *ldq, int *indxq, double *rho, int *cutpnt, double *z, double *dlamda, double *q2, int *ldq2, double *w, int *perm, int *givptr, int *givcol, double *givnum, int *indxp, int *indx, int *info);

F77_RET_I dlaed9_f77(int *k, int *kstart, int *kstop, int *n, double *d, double *q, int *ldq, double *rho, double *dlamda, double *w, double *s, int *lds, int *info);

F77_RET_I dlaeda_f77(int *n, int *tlvls, int *curlvl, int *curpbm, int *prmptr, int *perm, int *givptr, int *givcol, double *givnum, double *q, int *qptr, double *z, double *ztemp, int *info);

F77_RET_I dlaein_f77(logical *rightv, logical *noinit, int *n, double *h, int *ldh, double *wr, double *wi, double *vr, double *vi, double *b, int *ldb, double *work, double *eps3, double *smlnum, double *bignum, int *info);

F77_RET_I dlaev2_f77(double *a, double *b, double *c, double *rt1, double *rt2, double *cs1, double *sn1);

F77_RET_I dlaexc_f77(logical *wantq, int *n, double *t, int *ldt, double *q, int *ldq, int *j1, int *n1, int *n2, double *work, int *info);

F77_RET_I dlag2_f77(double *a, int *lda, double *b, int *ldb, double *safmin, double *scale1, double *scale2, double *wr1, double *wr2, double *wi);

F77_RET_I dlag2s_f77(int *m, int *n, double *a, int *lda, float *sa, int *ldsa, int *info);

F77_RET_I dlags2_f77(logical *upper, double *a1, double *a2, double *a3, double *b1, double *b2, double *b3, double *csu, double *snu, double *csv, double *snv, double *csq, double *snq);

F77_RET_I dlagtf_f77(int *n, double *a, double *lambda, double *b, double *c, double *tol, double *d, int *in, int *info);

F77_RET_I dlagtm_f77(const char *trans, int *n, int *nrhs, double *alpha, double *dl, double *d, double *du, double *x, int *ldx, double *beta, double *b, int *ldb);

F77_RET_I dlagts_f77(int *job, int *n, double *a, double *b, double *c, double *d, int *in, double *y, double *tol, int *info);

F77_RET_I dlagv2_f77(double *a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double *beta, double *csl, double *snl, double *csr, double *snr);

F77_RET_I dlahqr_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, double *h, int *ldh, double *wr, double *wi, int *iloz, int *ihiz, double *z, int *ldz, int *info);

F77_RET_I dlahr2_f77(int *n, int *k, int *nb, double *a, int *lda, double *tau, double *t, int *ldt, double *y, int *ldy);

F77_RET_I dlahrd_f77(int *n, int *k, int *nb, double *a, int *lda, double *tau, double *t, int *ldt, double *y, int *ldy);

F77_RET_I dlaic1_f77(int *job, int *j, double *x, double *sest, double *w, double *gamma, double *sestpr, double *s, double *c);

logical dlaisnan_f77(double *din1, double *din2);

F77_RET_I dlaln2_f77(logical *ltrans, int *na, int *nw, double *smin, double *ca, double *a, int *lda, double *d1, double *d2, double *b, int *ldb, double *wr, double *wi, double *x, int *ldx, double *scale, double *xnorm, int *info);

F77_RET_I dlals0_f77(int *icompq, int *nl, int *nr, int *sqre, int *nrhs, double *b, int *ldb, double *bx, int *ldbx, int *perm, int *givptr, int *givcol, int *ldgcol, double *givnum, int *ldgnum, double *poles, double *difl, double *difr, double *z, int *k, double *c, double *s, double *work, int *info);

F77_RET_I dlalsa_f77(int *icompq, int *smlsiz, int *n, int *nrhs, double *b, int *ldb, double *bx, int *ldbx, double *u, int *ldu, double *vt, int *k, double *difl, double *difr, double *z, double *poles, int *givptr, int *givcol, int *ldgcol, int *perm, double *givnum, double *c, double *s, double *work, int *iwork, int *info);

F77_RET_I dlalsd_f77(const char *uplo, int *smlsiz, int *n, int *nrhs, double *d, double *e, double *b, int *ldb, double *rcond, int *rank, double *work, int *iwork, int *info);

F77_RET_D dlamch_f77(const char *camch);

F77_RET_I dlamrg_f77(int *n1, int *n2, double *a, int *dtrd1, int *dtrd2, int *index);

F77_RET_I dlaneg_f77(int *n, double *d, double *lld, double *sigma, double *pivmin, int *r);

F77_RET_D dlangb_f77(const char *norm, int *n, int *kl, int *ku, double *ab, int *ldab, double *work);

F77_RET_D dlange_f77(const char *norm, int *m, int *n, double *a, int *lda, double *work);

F77_RET_D dlangt_f77(const char *norm, int *n, double *dl, double *d, double *du);

F77_RET_D dlanhs_f77(const char *norm, int *n, double *a, int *lda, double *work);

F77_RET_D dlansb_f77(const char *norm, const char *uplo, int *n, int *k, double *ab, int *ldab, double *work);

F77_RET_D dlansp_f77(const char *norm, const char *uplo, int *n, double *ap, double *work);

F77_RET_D dlanst_f77(const char *norm, int *n, double *d, double *e);

F77_RET_D dlansy_f77(const char *norm, const char *uplo, int *n, double *a, int *lda, double *work);

F77_RET_D dlantb_f77(const char *norm, const char *uplo, const char *diag, int *n, int *k, double *ab, int *ldab, double *work);

F77_RET_D dlantp_f77(const char *norm, const char *uplo, const char *diag, int *n, double *ap, double *work);

F77_RET_D dlantr_f77(const char *norm, const char *uplo, const char *diag, int *m, int *n, double *a, int *lda, double *work);

F77_RET_I dlanv2_f77(double *a, double *b, double *c, double *d, double *rt1r, double *rt1i, double *rt2r, double *rt2i, double *cs, double *sn);

F77_RET_I dlapll_f77(int *n, double *x, int *incx, double *y, int *incy, double *ssmin);

F77_RET_I dlapmt_f77(logical *forwrd, int *m, int *n, double *x, int *ldx, int *k);

F77_RET_D dlapy2_f77(double *x, double *y);

F77_RET_D dlapy3_f77(double *x, double *y, double *z);

F77_RET_I dlaqgb_f77(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, const char *equed);

F77_RET_I dlaqge_f77(int *m, int *n, double *a, int *lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, const char *equed);

F77_RET_I dlaqp2_f77(int *m, int *n, int *offset, double *a, int *lda, int *jpvt, double *tau, double *vn1, double *vn2, double *work);

F77_RET_I dlaqps_f77(int *m, int *n, int *offset, int *nb, int *kb, double *a, int *lda, int *jpvt, double *tau, double *vn1, double *vn2, double *auxv, double *f, int *ldf);

F77_RET_I dlaqr0_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, double *h, int *ldh, double *wr, double *wi, int *iloz, int *ihiz, double *z, int *ldz, double *work, int *lwork, int *info);

F77_RET_I dlaqr1_f77(int *n, double *h, int *ldh, double *sr1, double *si1, double *sr2, double *si2, double *v);

F77_RET_I dlaqr2_f77(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, double *h, int *ldh, int *iloz, int *ihiz, double *z, int *ldz, int *ns, int *nd, double *sr, double *si, double *v, int *ldv, int *nh, double *t, int *ldt, int *nv, double *wv, int *ldwv, double *work, int *lwork);

F77_RET_I dlaqr3_f77(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, double *h, int *ldh, int *iloz, int *ihiz, double *z, int *ldz, int *ns, int *nd, double *sr, double *si, double *v, int *ldv, int *nh, double *t, int *ldt, int *nv, double *wv, int *ldwv, double *work, int *lwork);

F77_RET_I dlaqr4_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, double *h, int *ldh, double *wr, double *wi, int *iloz, int *ihiz, double *z, int *ldz, double *work, int *lwork, int *info);

F77_RET_I dlaqr5_f77(logical *wantt, logical *wantz, int *kacc22, int *n, int *ktop, int *kbot, int *nshfts, double *sr, double *si, double *h, int *ldh, int *iloz, int *ihiz, double *z, int *ldz, double *v, int *ldv, double *u, int *ldu, int *nv, double *wv, int *ldwv, int *nh, double *wh, int *ldwh);

F77_RET_I dlaqsb_f77(const char *uplo, int *n, int *kd, double *ab, int *ldab, double *s, double *scond, double *amax, const char *equed);

F77_RET_I dlaqsp_f77(const char *uplo, int *n, double *ap, double *s, double *scond, double *amax, const char *equed);

F77_RET_I dlaqsy_f77(const char *uplo, int *n, double *a, int *lda, double *s, double *scond, double *amax, const char *equed);

F77_RET_I dlaqtr_f77(logical *ltran, logical *lfloat, int *n, double *t, int *ldt, double *b, double *w, double *scale, double *x, double *work, int *info);

F77_RET_I dlar1v_f77(int *n, int *b1, int *bn, double *lambda, double *d, double *l, double *ld, double *lld, double *pivmin, double *gaptol, double *z, logical *wantnc, int *negcnt, double *ztz, double *mingma, int *r, int *isuppz, double *nrminv, double *resid, double *rqcorr, double *work);

F77_RET_I dlar2v_f77(int *n, double *x, double *y, double *z, int *incx, double *c, double *s, int *incc);

F77_RET_I dlarf_f77(const char *side, int *m, int *n, double *v, int *incv, double *tau, double *c, int *ldc, double *work);

F77_RET_I dlarfb_f77(const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n, int *k, double *v, int *ldv, double *t, int *ldt, double *c, int *ldc, double *work, int *ldwork);

F77_RET_I dlarfg_f77(int *n, double *alpha, double *x, int *incx, double *tau);

F77_RET_I dlarft_f77(const char *direct, const char *storev, int *n, int *k, double *v, int *ldv, double *tau, double *t, int *ldt);

F77_RET_I dlarfx_f77(const char *side, int *m, int *n, double *v, double *tau, double *c, int *ldc, double *work);

F77_RET_I dlargv_f77(int *n, double *x, int *incx, double *y, int *incy, double *c, int *incc);

F77_RET_I dlarnv_f77(int *idist, int *iseed, int *n, double *x);

F77_RET_I dlarra_f77(int *n, double *d, double *e, double *e2, double *spltol, double *tnrm, int *nsplit, int *isplit, int *info);

F77_RET_I dlarrb_f77(int *n, double *d, double *lld, int *ifirst, int *ilast, double *rtol1, double *rtol2, int *offset, double *w, double *wgap, double *werr, double *work, int *iwork, double *pivmin, double *spdiam, int *twist, int *info);

F77_RET_I dlarrc_f77(const char *jobt, int *n, double *vl, double *vu, double *d, double *e, double *pivmin, int *eigcnt, int *lcnt, int *rcnt, int *info);

F77_RET_I dlarrd_f77(const char *range, const char *order, int *n, double *vl, double *vu, int *il, int *iu, double *gers, double *reltol, double *d, double *e, double *e2, double *pivmin, int *nsplit, int *isplit, int *m, double *w, double *werr, double *wl, double *wu, int *iblock, int *indexw, double *work, int *iwork, int *info);

F77_RET_I dlarre_f77(const char *range, int *n, double *vl, double *vu, int *il, int *iu, double *d, double *e, double *e2, double *rtol1, double *rtol2, double *spltol, int *nsplit, int *isplit, int *m, double *w, double *werr, double *wgap, int *iblock, int *indexw, double *gers, double *pivmin, double *work, int *iwork, int *info);

F77_RET_I dlarrf_f77(int *n, double *d, double *l, double *ld, int *clstrt, int *clend, double *w, double *wgap, double *werr, double *spdiam, double *clgapl, double *clgapr, double *pivmin, double *sigma, double *dplus, double *lplus, double *work, int *info);

F77_RET_I dlarrj_f77(int *n, double *d, double *e2, int *ifirst, int *ilast, double *rtol, int *offset, double *w, double *werr, double *work, int *iwork, double *pivmin, double *spdiam, int *info);

F77_RET_I dlarrk_f77(int *n, int *iw, double *gl, double *gu, double *d, double *e2, double *pivmin, double *reltol, double *w, double *werr, int *info);

F77_RET_I dlarrr_f77(int *n, double *d, double *e, int *info);

F77_RET_I dlarrv_f77(int *n, double *vl, double *vu, double *d, double *l, double *pivmin, int *isplit, int *m, int *dol, int *dou, double *minrgp, double *rtol1, double *rtol2, double *w, double *werr, double *wgap, int *iblock, int *indexw, double *gers, double *z, int *ldz, int *isuppz, double *work, int *iwork, int *info);

F77_RET_I dlartg_f77(double *f, double *g, double *cs, double *sn, double *r);

F77_RET_I dlartv_f77(int *n, double *x, int *incx, double *y, int *incy, double *c, double *s, int *incc);

F77_RET_I dlaruv_f77(int *iseed, int *n, double *x);

F77_RET_I dlarz_f77(const char *side, int *m, int *n, int *l, double *v, int *incv, double *tau, double *c, int *ldc, double *work);

F77_RET_I dlarzb_f77(const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n, int *k, int *l, double *v, int *ldv, double *t, int *ldt, double *c, int *ldc, double *work, int *ldwork);

F77_RET_I dlarzt_f77(const char *direct, const char *storev, int *n, int *k, double *v, int *ldv, double *tau, double *t, int *ldt);

F77_RET_I dlas2_f77(double *f, double *g, double *h, double *ssmin, double *ssmax);

F77_RET_I dlascl_f77(const char *type, int *kl, int *ku, double *cfrom, double *cto, int *m, int *n, double *a, int *lda, int *info);

F77_RET_I dlasd0_f77(int *n, int *sqre, double *d, double *e, double *u, int *ldu, double *vt, int *ldvt, int *smlsiz, int *iwork, double *work, int *info);

F77_RET_I dlasd1_f77(int *nl, int *nr, int *sqre, double *d, double *alpha, double *beta, double *u, int *ldu, double *vt, int *ldvt, int *idxq, int *iwork, double *work, int *info);

F77_RET_I dlasd2_f77(int *nl, int *nr, int *sqre, int *k, double *d, double *z, double *alpha, double *beta, double *u, int *ldu, double *vt, int *ldvt, double *dsigma, double *u2, int *ldu2, double *vt2, int *ldvt2, int *idxp, int *idx, int *idxc, int *idxq, int *coltyp, int *info);

F77_RET_I dlasd3_f77(int *nl, int *nr, int *sqre, int *k, double *d, double *q, int *ldq, double *dsigma, double *u, int *ldu, double *u2, int *ldu2, double *vt, int *ldvt, double *vt2, int *ldvt2, int *idxc, int *ctot, double *z, int *info);

F77_RET_I dlasd4_f77(int *n, int *i, double *d, double *z, double *delta, double *rho, double *sigma, double *work, int *info);

F77_RET_I dlasd5_f77(int *i, double *d, double *z, double *delta, double *rho, double *dsigma, double *work);

F77_RET_I dlasd6_f77(int *icompq, int *nl, int *nr, int *sqre, double *d, double *vf, double *vl, double *alpha, double *beta, int *idxq, int *perm, int *givptr, int *givcol, int *ldgcol, double *givnum, int *ldgnum, double *poles, double *difl, double *difr, double *z, int *k, double *c, double *s, double *work, int *iwork, int *info);

F77_RET_I dlasd7_f77(int *icompq, int *nl, int *nr, int *sqre, int *k, double *d, double *z, double *zw, double *vf, double *vfw, double *vl, double *vlw, double *alpha, double *beta, double *dsigma, int *idx, int *idxp, int *idxq, int *perm, int *givptr, int *givcol, int *ldgcol, double *givnum, int *ldgnum, double *c, double *s, int *info);

F77_RET_I dlasd8_f77(int *icompq, int *k, double *d, double *z, double *vf, double *vl, double *difl, double *difr, int *lddifr, double *dsigma, double *work, int *info);

F77_RET_I dlasda_f77(int *icompq, int *smlsiz, int *n, int *sqre, double *d, double *e, double *u, int *ldu, double *vt, int *k, double *difl, double *difr, double *z, double *poles, int *givptr, int *givcol, int *ldgcol, int *perm, double *givnum, double *c, double *s, double *work, int *iwork, int *info);

F77_RET_I dlasdq_f77(const char *uplo, int *sqre, int *n, int *ncvt, int *nru, int *ncc, double *d, double *e, double *vt, int *ldvt, double *u, int *ldu, double *c, int *ldc, double *work, int *info);

F77_RET_I dlasdt_f77(int *n, int *lvl, int *nd, int *inode, int *ndiml, int *ndimr, int *msub);

F77_RET_I dlaset_f77(const char *uplo, int *m, int *n, double *alpha, double *beta, double *a, int *lda);

F77_RET_I dlasq1_f77(int *n, double *d, double *e, double *work, int *info);

F77_RET_I dlasq2_f77(int *n, double *z, int *info);

F77_RET_I dlasq3_f77(int *i0, int *n0, double *z, int *pp, double *dmin, double *sigma, double *desig, double *qmax, int *nfail, int *iter, int *ndiv, logical *ieee);

F77_RET_I dlasq4_f77(int *i0, int *n0, double *z, int *pp, int *n0in, double *dmin, double *dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *tau, int *ttype);

F77_RET_I dlasq5_f77(int *i0, int *n0, double *z, int *pp, double *tau, double *dmin, double *dmin1, double *dmin2, double *dn, double *dnm1, double *dnm2, logical *ieee);

F77_RET_I dlasq6_f77(int *i0, int *n0, double *z, int *pp, double *dmin, double *dmin1, double *dmin2, double *dn, double *dnm1, double *dnm2);

F77_RET_I dlasr_f77(const char *side, const char *pivot, const char *direct, int *m, int *n, double *c, double *s, double *a, int *lda);

F77_RET_I dlasrt_f77(const char *id, int *n, double *d, int *info);

F77_RET_I dlassq_f77(int *n, double *x, int *incx, double *scale, double *sumsq);

F77_RET_I dlasv2_f77(double *f, double *g, double *h, double *ssmin, double *ssmax, double *snr, double *csr, double *snl, double *csl);

F77_RET_I dlaswp_f77(int *n, double *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);

F77_RET_I dlasy2_f77(logical *ltranl, logical *ltranr, int *isgn, int *n1, int *n2, double *tl, int *ldtl, double *tr, int *ldtr, double *b, int *ldb, double *scale, double *x, int *ldx, double *xnorm, int *info);

F77_RET_I dlasyf_f77(const char *uplo, int *n, int *nb, int *kb, double *a, int *lda, int *ipiv, double *w, int *ldw, int *info);

F77_RET_I dlatbs_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, int *kd, double *ab, int *ldab, double *x, double *scale, double *cnorm, int *info);

F77_RET_I dlatdf_f77(int *ijob, int *n, double *z, int *ldz, double *rhs, double *rdsum, double *rdscal, int *ipiv, int *jpiv);

F77_RET_I dlatps_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, double *ap, double *x, double *scale, double *cnorm, int *info);

F77_RET_I dlatrd_f77(const char *uplo, int *n, int *nb, double *a, int *lda, double *e, double *tau, double *w, int *ldw);

F77_RET_I dlatrs_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, double *a, int *lda, double *x, double *scale, double *cnorm, int *info);

F77_RET_I dlatrz_f77(int *m, int *n, int *l, double *a, int *lda, double *tau, double *work);

F77_RET_I dlatzm_f77(const char *side, int *m, int *n, double *v, int *incv, double *tau, double *c1, double *c2, int *ldc, double *work);

F77_RET_I dlauu2_f77(const char *uplo, int *n, double *a, int *lda, int *info);

F77_RET_I dlauum_f77(const char *uplo, int *n, double *a, int *lda, int *info);

F77_RET_I dlazq3_f77(int *i0, int *n0, double *z, int *pp, double *dmin, double *sigma, double *desig, double *qmax, int *nfail, int *iter, int *ndiv, logical *ieee, int *ttype, double *dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *tau);

F77_RET_I dlazq4_f77(int *i0, int *n0, double *z, int *pp, int *n0in, double *dmin, double *dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *tau, int *ttype, double *g);

F77_RET_I dopgtr_f77(const char *uplo, int *n, double *ap, double *tau, double *q, int *ldq, double *work, int *info);

F77_RET_I dopmtr_f77(const char *side, const char *uplo, const char *trans, int *m, int *n, double *ap, double *tau, double *c, int *ldc, double *work, int *info);

F77_RET_I dorg2l_f77(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *info);

F77_RET_I dorg2r_f77(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *info);

F77_RET_I dorgbr_f77(const char *vect, int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dorghr_f77(int *n, int *ilo, int *ihi, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dorgl2_f77(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *info);

F77_RET_I dorglq_f77(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dorgql_f77(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dorgqr_f77(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dorgr2_f77(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *info);

F77_RET_I dorgrq_f77(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dorgtr_f77(const char *uplo, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_I dorm2l_f77(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *info);

F77_RET_I dorm2r_f77(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *info);

F77_RET_I dormbr_f77(const char *vect, const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

F77_RET_I dormhr_f77(const char *side, const char *trans, int *m, int *n, int *ilo, int *ihi, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

F77_RET_I dorml2_f77(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *info);

F77_RET_I dormlq_f77(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

F77_RET_I dormql_f77(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

F77_RET_I dormqr_f77(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

F77_RET_I dormr2_f77(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *info);

F77_RET_I dormr3_f77(const char *side, const char *trans, int *m, int *n, int *k, int *l, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *info);

F77_RET_I dormrq_f77(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

F77_RET_I dormrz_f77(const char *side, const char *trans, int *m, int *n, int *k, int *l, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

F77_RET_I dormtr_f77(const char *side, const char *uplo, const char *trans, int *m, int *n, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

F77_RET_I dpbcon_f77(const char *uplo, int *n, int *kd, double *ab, int *ldab, double *anorm, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dpbequ_f77(const char *uplo, int *n, int *kd, double *ab, int *ldab, double *s, double *scond, double *amax, int *info);

F77_RET_I dpbrfs_f77(const char *uplo, int *n, int *kd, int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dpbstf_f77(const char *uplo, int *n, int *kd, double *ab, int *ldab, int *info);

F77_RET_I dpbsv_f77(const char *uplo, int *n, int *kd, int *nrhs, double *ab, int *ldab, double *b, int *ldb, int *info);

F77_RET_I dpbsvx_f77(const char *fact, const char *uplo, int *n, int *kd, int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, const char *equed, double *s, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dpbtf2_f77(const char *uplo, int *n, int *kd, double *ab, int *ldab, int *info);

F77_RET_I dpbtrf_f77(const char *uplo, int *n, int *kd, double *ab, int *ldab, int *info);

F77_RET_I dpbtrs_f77(const char *uplo, int *n, int *kd, int *nrhs, double *ab, int *ldab, double *b, int *ldb, int *info);

F77_RET_I dpocon_f77(const char *uplo, int *n, double *a, int *lda, double *anorm, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dpoequ_f77(int *n, double *a, int *lda, double *s, double *scond, double *amax, int *info);

F77_RET_I dporfs_f77(const char *uplo, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dposv_f77(const char *uplo, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);

F77_RET_I dposvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, const char *equed, double *s, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dpotf2_f77(const char *uplo, int *n, double *a, int *lda, int *info);

F77_RET_I dpotrf_f77(const char *uplo, int *n, double *a, int *lda, int *info);

F77_RET_I dpotri_f77(const char *uplo, int *n, double *a, int *lda, int *info);

F77_RET_I dpotrs_f77(const char *uplo, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);

F77_RET_I dppcon_f77(const char *uplo, int *n, double *ap, double *anorm, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dppequ_f77(const char *uplo, int *n, double *ap, double *s, double *scond, double *amax, int *info);

F77_RET_I dpprfs_f77(const char *uplo, int *n, int *nrhs, double *ap, double *afp, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dppsv_f77(const char *uplo, int *n, int *nrhs, double *ap, double *b, int *ldb, int *info);

F77_RET_I dppsvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, double *ap, double *afp, const char *equed, double *s, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dpptrf_f77(const char *uplo, int *n, double *ap, int *info);

F77_RET_I dpptri_f77(const char *uplo, int *n, double *ap, int *info);

F77_RET_I dpptrs_f77(const char *uplo, int *n, int *nrhs, double *ap, double *b, int *ldb, int *info);

F77_RET_I dptcon_f77(int *n, double *d, double *e, double *anorm, double *rcond, double *work, int *info);

F77_RET_I dpteqr_f77(const char *compz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *info);

F77_RET_I dptrfs_f77(int *n, int *nrhs, double *d, double *e, double *df, double *ef, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *info);

F77_RET_I dptsv_f77(int *n, int *nrhs, double *d, double *e, double *b, int *ldb, int *info);

F77_RET_I dptsvx_f77(const char *fact, int *n, int *nrhs, double *d, double *e, double *df, double *ef, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *info);

F77_RET_I dpttrf_f77(int *n, double *d, double *e, int *info);

F77_RET_I dpttrs_f77(int *n, int *nrhs, double *d, double *e, double *b, int *ldb, int *info);

F77_RET_I dptts2_f77(int *n, int *nrhs, double *d, double *e, double *b, int *ldb);

F77_RET_I drscl_f77(int *n, double *sa, double *sx, int *incx);

F77_RET_I dsbev_f77(const char *jobz, const char *uplo, int *n, int *kd, double *ab, int *ldab, double *w, double *z, int *ldz, double *work, int *info);

F77_RET_I dsbevd_f77(const char *jobz, const char *uplo, int *n, int *kd, double *ab, int *ldab, double *w, double *z, int *ldz, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dsbevx_f77(const char *jobz, const char *range, const char *uplo, int *n, int *kd, double *ab, int *ldab, double *q, int *ldq, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, double *work, int *iwork, int *ifail, int *info);

F77_RET_I dsbgst_f77(const char *vect, const char *uplo, int *n, int *ka, int *kb, double *ab, int *ldab, double *bb, int *ldbb, double *x, int *ldx, double *work, int *info);

F77_RET_I dsbgv_f77(const char *jobz, const char *uplo, int *n, int *ka, int *kb, double *ab, int *ldab, double *bb, int *ldbb, double *w, double *z, int *ldz, double *work, int *info);

F77_RET_I dsbgvd_f77(const char *jobz, const char *uplo, int *n, int *ka, int *kb, double *ab, int *ldab, double *bb, int *ldbb, double *w, double *z, int *ldz, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dsbgvx_f77(const char *jobz, const char *range, const char *uplo, int *n, int *ka, int *kb, double *ab, int *ldab, double *bb, int *ldbb, double *q, int *ldq, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, double *work, int *iwork, int *ifail, int *info);

F77_RET_I dsbtrd_f77(const char *vect, const char *uplo, int *n, int *kd, double *ab, int *ldab, double *d, double *e, double *q, int *ldq, double *work, int *info);

F77_RET_I dsgesv_f77(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, double *x, int *ldx, double *work, float *swork, int *iter, int *info);

F77_RET_I dspcon_f77(const char *uplo, int *n, double *ap, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dspev_f77(const char *jobz, const char *uplo, int *n, double *ap, double *w, double *z, int *ldz, double *work, int *info);

F77_RET_I dspevd_f77(const char *jobz, const char *uplo, int *n, double *ap, double *w, double *z, int *ldz, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dspevx_f77(const char *jobz, const char *range, const char *uplo, int *n, double *ap, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, double *work, int *iwork, int *ifail, int *info);

F77_RET_I dspgst_f77(int *itype, const char *uplo, int *n, double *ap, double *bp, int *info);

F77_RET_I dspgv_f77(int *itype, const char *jobz, const char *uplo, int *n, double *ap, double *bp, double *w, double *z, int *ldz, double *work, int *info);

F77_RET_I dspgvd_f77(int *itype, const char *jobz, const char *uplo, int *n, double *ap, double *bp, double *w, double *z, int *ldz, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dspgvx_f77(int *itype, const char *jobz, const char *range, const char *uplo, int *n, double *ap, double *bp, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, double *work, int *iwork, int *ifail, int *info);

F77_RET_I dsprfs_f77(const char *uplo, int *n, int *nrhs, double *ap, double *afp, int *ipiv, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dspsv_f77(const char *uplo, int *n, int *nrhs, double *ap, int *ipiv, double *b, int *ldb, int *info);

F77_RET_I dspsvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, double *ap, double *afp, int *ipiv, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dsptrd_f77(const char *uplo, int *n, double *ap, double *d, double *e, double *tau, int *info);

F77_RET_I dsptrf_f77(const char *uplo, int *n, double *ap, int *ipiv, int *info);

F77_RET_I dsptri_f77(const char *uplo, int *n, double *ap, int *ipiv, double *work, int *info);

F77_RET_I dsptrs_f77(const char *uplo, int *n, int *nrhs, double *ap, int *ipiv, double *b, int *ldb, int *info);

F77_RET_I dstebz_f77(const char *range, const char *order, int *n, double *vl, double *vu, int *il, int *iu, double *abstol, double *d, double *e, int *m, int *nsplit, double *w, int *iblock, int *isplit, double *work, int *iwork, int *info);

F77_RET_I dstedc_f77(const char *compz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dstegr_f77(const char *jobz, const char *range, int *n, double *d, double *e, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dstein_f77(int *n, double *d, double *e, int *m, double *w, int *iblock, int *isplit, double *z, int *ldz, double *work, int *iwork, int *ifail, int *info);

F77_RET_I dstemr_f77(const char *jobz, const char *range, int *n, double *d, double *e, double *vl, double *vu, int *il, int *iu, int *m, double *w, double *z, int *ldz, int *nzc, int *isuppz, logical *tryrac, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dsteqr_f77(const char *compz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *info);

F77_RET_I dsterf_f77(int *n, double *d, double *e, int *info);

F77_RET_I dstev_f77(const char *jobz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *info);

F77_RET_I dstevd_f77(const char *jobz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dstevr_f77(const char *jobz, const char *range, int *n, double *d, double *e, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dstevx_f77(const char *jobz, const char *range, int *n, double *d, double *e, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, double *work, int *iwork, int *ifail, int *info);

F77_RET_I dsycon_f77(const char *uplo, int *n, double *a, int *lda, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dsyev_f77(const char *jobz, const char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);

F77_RET_I dsyevd_f77(const char *jobz, const char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dsyevr_f77(const char *jobz, const char *range, const char *uplo, int *n, double *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dsyevx_f77(const char *jobz, const char *range, const char *uplo, int *n, double *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, double *work, int *lwork, int *iwork, int *ifail, int *info);

F77_RET_I dsygs2_f77(int *itype, const char *uplo, int *n, double *a, int *lda, double *b, int *ldb, int *info);

F77_RET_I dsygst_f77(int *itype, const char *uplo, int *n, double *a, int *lda, double *b, int *ldb, int *info);

F77_RET_I dsygv_f77(int *itype, const char *jobz, const char *uplo, int *n, double *a, int *lda, double *b, int *ldb, double *w, double *work, int *lwork, int *info);

F77_RET_I dsygvd_f77(int *itype, const char *jobz, const char *uplo, int *n, double *a, int *lda, double *b, int *ldb, double *w, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dsygvx_f77(int *itype, const char *jobz, const char *range, const char *uplo, int *n, double *a, int *lda, double *b, int *ldb, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, double *work, int *lwork, int *iwork, int *ifail, int *info);

F77_RET_I dsyrfs_f77(const char *uplo, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dsysv_f77(const char *uplo, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, double *work, int *lwork, int *info);

F77_RET_I dsysvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *lwork, int *iwork, int *info);

F77_RET_I dsytd2_f77(const char *uplo, int *n, double *a, int *lda, double *d, double *e, double *tau, int *info);

F77_RET_I dsytf2_f77(const char *uplo, int *n, double *a, int *lda, int *ipiv, int *info);

F77_RET_I dsytrd_f77(const char *uplo, int *n, double *a, int *lda, double *d, double *e, double *tau, double *work, int *lwork, int *info);

F77_RET_I dsytrf_f77(const char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

F77_RET_I dsytri_f77(const char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *info);

F77_RET_I dsytrs_f77(const char *uplo, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

F77_RET_I dtbcon_f77(const char *norm, const char *uplo, const char *diag, int *n, int *kd, double *ab, int *ldab, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dtbrfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *kd, int *nrhs, double *ab, int *ldab, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dtbtrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *kd, int *nrhs, double *ab, int *ldab, double *b, int *ldb, int *info);

F77_RET_I dtgevc_f77(const char *side, const char *howmny, logical *select, int *n, double *s, int *lds, double *p, int *ldp, double *vl, int *ldvl, double *vr, int *ldvr, int *mm, int *m, double *work, int *info);

F77_RET_I dtgex2_f77(logical *wantq, logical *wantz, int *n, double *a, int *lda, double *b, int *ldb, double *q, int *ldq, double *z, int *ldz, int *j1, int *n1, int *n2, double *work, int *lwork, int *info);

F77_RET_I dtgexc_f77(logical *wantq, logical *wantz, int *n, double *a, int *lda, double *b, int *ldb, double *q, int *ldq, double *z, int *ldz, int *ifst, int *ilst, double *work, int *lwork, int *info);

F77_RET_I dtgsen_f77(int *ijob, logical *wantq, logical *wantz, logical *select, int *n, double *a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double *beta, double *q, int *ldq, double *z, int *ldz, int *m, double *pl, double *pr, double *dif, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dtgsja_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *p, int *n, int *k, int *l, double *a, int *lda, double *b, int *ldb, double *tola, double *tolb, double *alpha, double *beta, double *u, int *ldu, double *v, int *ldv, double *q, int *ldq, double *work, int *ncycle, int *info);

F77_RET_I dtgsna_f77(const char *job, const char *howmny, logical *select, int *n, double *a, int *lda, double *b, int *ldb, double *vl, int *ldvl, double *vr, int *ldvr, double *s, double *dif, int *mm, int *m, double *work, int *lwork, int *iwork, int *info);

F77_RET_I dtgsy2_f77(const char *trans, int *ijob, int *m, int *n, double *a, int *lda, double *b, int *ldb, double *c, int *ldc, double *d, int *ldd, double *e, int *lde, double *f, int *ldf, double *scale, double *rdsum, double *rdscal, int *iwork, int *pq, int *info);

F77_RET_I dtgsyl_f77(const char *trans, int *ijob, int *m, int *n, double *a, int *lda, double *b, int *ldb, double *c, int *ldc, double *d, int *ldd, double *e, int *lde, double *f, int *ldf, double *scale, double *dif, double *work, int *lwork, int *iwork, int *info);

F77_RET_I dtpcon_f77(const char *norm, const char *uplo, const char *diag, int *n, double *ap, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dtprfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, double *ap, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dtptri_f77(const char *uplo, const char *diag, int *n, double *ap, int *info);

F77_RET_I dtptrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, double *ap, double *b, int *ldb, int *info);

F77_RET_I dtrcon_f77(const char *norm, const char *uplo, const char *diag, int *n, double *a, int *lda, double *rcond, double *work, int *iwork, int *info);

F77_RET_I dtrevc_f77(const char *side, const char *howmny, logical *select, int *n, double *t, int *ldt, double *vl, int *ldvl, double *vr, int *ldvr, int *mm, int *m, double *work, int *info);

F77_RET_I dtrexc_f77(const char *compq, int *n, double *t, int *ldt, double *q, int *ldq, int *ifst, int *ilst, double *work, int *info);

F77_RET_I dtrrfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);

F77_RET_I dtrsen_f77(const char *job, const char *compq, logical *select, int *n, double *t, int *ldt, double *q, int *ldq, double *wr, double *wi, int *m, double *s, double *sep, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I dtrsna_f77(const char *job, const char *howmny, logical *select, int *n, double *t, int *ldt, double *vl, int *ldvl, double *vr, int *ldvr, double *s, double *sep, int *mm, int *m, double *work, int *ldwork, int *iwork, int *info);

F77_RET_I dtrsyl_f77(const char *trana, const char *tranb, int *isgn, int *m, int *n, double *a, int *lda, double *b, int *ldb, double *c, int *ldc, double *scale, int *info);

F77_RET_I dtrti2_f77(const char *uplo, const char *diag, int *n, double *a, int *lda, int *info);

F77_RET_I dtrtri_f77(const char *uplo, const char *diag, int *n, double *a, int *lda, int *info);

F77_RET_I dtrtrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);

F77_RET_I dtzrqf_f77(int *m, int *n, double *a, int *lda, double *tau, int *info);

F77_RET_I dtzrzf_f77(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

F77_RET_D dzsum1_f77(int *n, doublecomplex *cx, int *incx);

F77_RET_I icmax1_f77(int *n, singlecomplex *cx, int *incx);

F77_RET_I ieeeck_f77(int *ispec, float *zero, float *one);

F77_RET_I ilaenv_f77(int *ispec, const char *name, const char *opts, int *n1, int *n2, int *n3, int *n4);

F77_RET_I ilaver_f77(int *vers_major, int *vers_minor, int *vers_patch);

F77_RET_I iparmq_f77(int *ispec, const char *name, const char *opts, int *n, int *ilo, int *ihi, int *lwork);

F77_RET_I izmax1_f77(int *n, doublecomplex *cx, int *incx);

logical lsamen_f77(int *n, const char *ca, const char *cb);

F77_RET_I sbdsdc_f77(const char *uplo, const char *compq, int *n, float *d, float *e, float *u, int *ldu, float *vt, int *ldvt, float *q, int *iq, float *work, int *iwork, int *info);

F77_RET_I sbdsqr_f77(const char *uplo, int *n, int *ncvt, int *nru, int *ncc, float *d, float *e, float *vt, int *ldvt, float *u, int *ldu, float *c, int *ldc, float *work, int *info);

F77_RET_F scsum1_f77(int *n, singlecomplex *cx, int *incx);

F77_RET_I sdisna_f77(const char *job, int *m, int *n, float *d, float *sep, int *info);

F77_RET_I sgbbrd_f77(const char *vect, int *m, int *n, int *ncc, int *kl, int *ku, float *ab, int *ldab, float *d, float *e, float *q, int *ldq, float *pt, int *ldpt, float *c, int *ldc, float *work, int *info);

F77_RET_I sgbcon_f77(const char *norm, int *n, int *kl, int *ku, float *ab, int *ldab, int *ipiv, float *anorm, float *rcond, float *work, int *iwork, int *info);

F77_RET_I sgbequ_f77(int *m, int *n, int *kl, int *ku, float *ab, int *ldab, float *r, float *c, float *rowcnd, float *colcnd, float *amax, int *info);

F77_RET_I sgbrfs_f77(const char *trans, int *n, int *kl, int *ku, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I sgbsv_f77(int *n, int *kl, int *ku, int *nrhs, float *ab, int *ldab, int *ipiv, float *b, int *ldb, int *info);

F77_RET_I sgbsvx_f77(const char *fact, const char *trans, int *n, int *kl, int *ku, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, int *ipiv, const char *equed, float *r, float *c, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I sgbtf2_f77(int *m, int *n, int *kl, int *ku, float *ab, int *ldab, int *ipiv, int *info);

F77_RET_I sgbtrf_f77(int *m, int *n, int *kl, int *ku, float *ab, int *ldab, int *ipiv, int *info);

F77_RET_I sgbtrs_f77(const char *trans, int *n, int *kl, int *ku, int *nrhs, float *ab, int *ldab, int *ipiv, float *b, int *ldb, int *info);

F77_RET_I sgebak_f77(const char *job, const char *side, int *n, int *ilo, int *ihi, float *scale, int *m, float *v, int *ldv, int *info);

F77_RET_I sgebal_f77(const char *job, int *n, float *a, int *lda, int *ilo, int *ihi, float *scale, int *info);

F77_RET_I sgebd2_f77(int *m, int *n, float *a, int *lda, float *d, float *e, float *tauq, float *taup, float *work, int *info);

F77_RET_I sgebrd_f77(int *m, int *n, float *a, int *lda, float *d, float *e, float *tauq, float *taup, float *work, int *lwork, int *info);

F77_RET_I sgecon_f77(const char *norm, int *n, float *a, int *lda, float *anorm, float *rcond, float *work, int *iwork, int *info);

F77_RET_I sgeequ_f77(int *m, int *n, float *a, int *lda, float *r, float *c, float *rowcnd, float *colcnd, float *amax, int *info);

F77_RET_I sgees_f77(const char *jobvs, const char *sort, L_fp select, int *n, float *a, int *lda, int *sdim, float *wr, float *wi, float *vs, int *ldvs, float *work, int *lwork, logical *bwork, int *info);

F77_RET_I sgeesx_f77(const char *jobvs, const char *sort, L_fp select, const char *sense, int *n, float *a, int *lda, int *sdim, float *wr, float *wi, float *vs, int *ldvs, float *rconde, float *rcondv, float *work, int *lwork, int *iwork, int *liwork, logical *bwork, int *info);

F77_RET_I sgeev_f77(const char *jobvl, const char *jobvr, int *n, float *a, int *lda, float *wr, float *wi, float *vl, int *ldvl, float *vr, int *ldvr, float *work, int *lwork, int *info);

F77_RET_I sgeevx_f77(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, int *n, float *a, int *lda, float *wr, float *wi, float *vl, int *ldvl, float *vr, int *ldvr, int *ilo, int *ihi, float *scale, float *abnrm, float *rconde, float *rcondv, float *work, int *lwork, int *iwork, int *info);

F77_RET_I sgegs_f77(const char *jobvsl, const char *jobvsr, int *n, float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *work, int *lwork, int *info);

F77_RET_I sgegv_f77(const char *jobvl, const char *jobvr, int *n, float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *vl, int *ldvl, float *vr, int *ldvr, float *work, int *lwork, int *info);

F77_RET_I sgehd2_f77(int *n, int *ilo, int *ihi, float *a, int *lda, float *tau, float *work, int *info);

F77_RET_I sgehrd_f77(int *n, int *ilo, int *ihi, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sgelq2_f77(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);

F77_RET_I sgelqf_f77(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sgels_f77(const char *trans, int *m, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, float *work, int *lwork, int *info);

F77_RET_I sgelsd_f77(int *m, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, float *s, float *rcond, int *rank, float *work, int *lwork, int *iwork, int *info);

F77_RET_I sgelss_f77(int *m, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, float *s, float *rcond, int *rank, float *work, int *lwork, int *info);

F77_RET_I sgelsx_f77(int *m, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, int *jpvt, float *rcond, int *rank, float *work, int *info);

F77_RET_I sgelsy_f77(int *m, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, int *jpvt, float *rcond, int *rank, float *work, int *lwork, int *info);

F77_RET_I sgeql2_f77(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);

F77_RET_I sgeqlf_f77(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sgeqp3_f77(int *m, int *n, float *a, int *lda, int *jpvt, float *tau, float *work, int *lwork, int *info);

F77_RET_I sgeqpf_f77(int *m, int *n, float *a, int *lda, int *jpvt, float *tau, float *work, int *info);

F77_RET_I sgeqr2_f77(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);

F77_RET_I sgeqrf_f77(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sgerfs_f77(const char *trans, int *n, int *nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I sgerq2_f77(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);

F77_RET_I sgerqf_f77(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sgesc2_f77(int *n, float *a, int *lda, float *rhs, int *ipiv, int *jpiv, float *scale);

F77_RET_I sgesdd_f77(const char *jobz, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *iwork, int *info);

F77_RET_I sgesv_f77(int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);

F77_RET_I sgesvd_f77(const char *jobu, const char *jobvt, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *info);

F77_RET_I sgesvx_f77(const char *fact, const char *trans, int *n, int *nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, const char *equed, float *r, float *c, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I sgetc2_f77(int *n, float *a, int *lda, int *ipiv, int *jpiv, int *info);

F77_RET_I sgetf2_f77(int *m, int *n, float *a, int *lda, int *ipiv, int *info);

F77_RET_I sgetrf_f77(int *m, int *n, float *a, int *lda, int *ipiv, int *info);

F77_RET_I sgetri_f77(int *n, float *a, int *lda, int *ipiv, float *work, int *lwork, int *info);

F77_RET_I sgetrs_f77(const char *trans, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);

F77_RET_I sggbak_f77(const char *job, const char *side, int *n, int *ilo, int *ihi, float *lscale, float *rscale, int *m, float *v, int *ldv, int *info);

F77_RET_I sggbal_f77(const char *job, int *n, float *a, int *lda, float *b, int *ldb, int *ilo, int *ihi, float *lscale, float *rscale, float *work, int *info);

F77_RET_I sgges_f77(const char *jobvsl, const char *jobvsr, const char *sort, L_fp selctg, int *n, float *a, int *lda, float *b, int *ldb, int *sdim, float *alphar, float *alphai, float *beta, float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *work, int *lwork, logical *bwork, int *info);

F77_RET_I sggesx_f77(const char *jobvsl, const char *jobvsr, const char *sort, L_fp selctg, const char *sense, int *n, float *a, int *lda, float *b, int *ldb, int *sdim, float *alphar, float *alphai, float *beta, float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *rconde, float *rcondv, float *work, int *lwork, int *iwork, int *liwork, logical *bwork, int *info);

F77_RET_I sggev_f77(const char *jobvl, const char *jobvr, int *n, float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *vl, int *ldvl, float *vr, int *ldvr, float *work, int *lwork, int *info);

F77_RET_I sggevx_f77(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, int *n, float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *vl, int *ldvl, float *vr, int *ldvr, int *ilo, int *ihi, float *lscale, float *rscale, float *abnrm, float *bbnrm, float *rconde, float *rcondv, float *work, int *lwork, int *iwork, logical *bwork, int *info);

F77_RET_I sggglm_f77(int *n, int *m, int *p, float *a, int *lda, float *b, int *ldb, float *d, float *x, float *y, float *work, int *lwork, int *info);

F77_RET_I sgghrd_f77(const char *compq, const char *compz, int *n, int *ilo, int *ihi, float *a, int *lda, float *b, int *ldb, float *q, int *ldq, float *z, int *ldz, int *info);

F77_RET_I sgglse_f77(int *m, int *n, int *p, float *a, int *lda, float *b, int *ldb, float *c, float *d, float *x, float *work, int *lwork, int *info);

F77_RET_I sggqrf_f77(int *n, int *m, int *p, float *a, int *lda, float *taua, float *b, int *ldb, float *taub, float *work, int *lwork, int *info);

F77_RET_I sggrqf_f77(int *m, int *p, int *n, float *a, int *lda, float *taua, float *b, int *ldb, float *taub, float *work, int *lwork, int *info);

F77_RET_I sggsvd_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *n, int *p, int *k, int *l, float *a, int *lda, float *b, int *ldb, float *alpha, float *beta, float *u, int *ldu, float *v, int *ldv, float *q, int *ldq, float *work, int *iwork, int *info);

F77_RET_I sggsvp_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *p, int *n, float *a, int *lda, float *b, int *ldb, float *tola, float *tolb, int *k, int *l, float *u, int *ldu, float *v, int *ldv, float *q, int *ldq, int *iwork, float *tau, float *work, int *info);

F77_RET_I sgtcon_f77(const char *norm, int *n, float *dl, float *d, float *du, float *du2, int *ipiv, float *anorm, float *rcond, float *work, int *iwork, int *info);

F77_RET_I sgtrfs_f77(const char *trans, int *n, int *nrhs, float *dl, float *d, float *du, float *dlf, float *df, float *duf, float *du2, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I sgtsv_f77(int *n, int *nrhs, float *dl, float *d, float *du, float *b, int *ldb, int *info);

F77_RET_I sgtsvx_f77(const char *fact, const char *trans, int *n, int *nrhs, float *dl, float *d, float *du, float *dlf, float *df, float *duf, float *du2, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I sgttrf_f77(int *n, float *dl, float *d, float *du, float *du2, int *ipiv, int *info);

F77_RET_I sgttrs_f77(const char *trans, int *n, int *nrhs, float *dl, float *d, float *du, float *du2, int *ipiv, float *b, int *ldb, int *info);

F77_RET_I sgtts2_f77(int *itrans, int *n, int *nrhs, float *dl, float *d, float *du, float *du2, int *ipiv, float *b, int *ldb);

F77_RET_I shgeqz_f77(const char *job, const char *compq, const char *compz, int *n, int *ilo, int *ihi, float *h, int *ldh, float *t, int *ldt, float *alphar, float *alphai, float *beta, float *q, int *ldq, float *z, int *ldz, float *work, int *lwork, int *info);

F77_RET_I shsein_f77(const char *side, const char *eigsrc, const char *initv, logical *select, int *n, float *h, int *ldh, float *wr, float *wi, float *vl, int *ldvl, float *vr, int *ldvr, int *mm, int *m, float *work, int *ifaill, int *ifailr, int *info);

F77_RET_I shseqr_f77(const char *job, const char *compz, int *n, int *ilo, int *ihi, float *h, int *ldh, float *wr, float *wi, float *z, int *ldz, float *work, int *lwork, int *info);

logical sisnan_f77(float *sin);

F77_RET_I slabad_f77(float *small, float *large);

F77_RET_I slabrd_f77(int *m, int *n, int *nb, float *a, int *lda, float *d, float *e, float *tauq, float *taup, float *x, int *ldx, float *y, int *ldy);

F77_RET_I slacn2_f77(int *n, float *v, float *x, int *isgn, float *est, int *kase, int *isave);

F77_RET_I slacon_f77(int *n, float *v, float *x, int *isgn, float *est, int *kase);

F77_RET_I slacpy_f77(const char *uplo, int *m, int *n, float *a, int *lda, float *b, int *ldb);

F77_RET_I sladiv_f77(float *a, float *b, float *c, float *d, float *p, float *q);

F77_RET_I slae2_f77(float *a, float *b, float *c, float *rt1, float *rt2);

F77_RET_I slaebz_f77(int *ijob, int *nitmax, int *n, int *mmax, int *minp, int *nbmin, float *abstol, float *reltol, float *pivmin, float *d, float *e, float *e2, int *nval, float *ab, float *c, int *mout, int *nab, float *work, int *iwork, int *info);

F77_RET_I slaed0_f77(int *icompq, int *qsiz, int *n, float *d, float *e, float *q, int *ldq, float *qstore, int *ldqs, float *work, int *iwork, int *info);

F77_RET_I slaed1_f77(int *n, float *d, float *q, int *ldq, int *indxq, float *rho, int *cutpnt, float *work, int *iwork, int *info);

F77_RET_I slaed2_f77(int *k, int *n, int *n1, float *d, float *q, int *ldq, int *indxq, float *rho, float *z, float *dlamda, float *w, float *q2, int *indx, int *indxc, int *indxp, int *coltyp, int *info);

F77_RET_I slaed3_f77(int *k, int *n, int *n1, float *d, float *q, int *ldq, float *rho, float *dlamda, float *q2, int *indx, int *ctot, float *w, float *s, int *info);

F77_RET_I slaed4_f77(int *n, int *i, float *d, float *z, float *delta, float *rho, float *dlam, int *info);

F77_RET_I slaed5_f77(int *i, float *d, float *z, float *delta, float *rho, float *dlam);

F77_RET_I slaed6_f77(int *kniter, logical *orgati, float *rho, float *d, float *z, float *finit, float *tau, int *info);

F77_RET_I slaed7_f77(int *icompq, int *n, int *qsiz, int *tlvls, int *curlvl, int *curpbm, float *d, float *q, int *ldq, int *indxq, float *rho, int *cutpnt, float *qstore, int *qptr, int *prmptr, int *perm, int *givptr, int *givcol, float *givnum, float *work, int *iwork, int *info);

F77_RET_I slaed8_f77(int *icompq, int *k, int *n, int *qsiz, float *d, float *q, int *ldq, int *indxq, float *rho, int *cutpnt, float *z, float *dlamda, float *q2, int *ldq2, float *w, int *perm, int *givptr, int *givcol, float *givnum, int *indxp, int *indx, int *info);

F77_RET_I slaed9_f77(int *k, int *kstart, int *kstop, int *n, float *d, float *q, int *ldq, float *rho, float *dlamda, float *w, float *s, int *lds, int *info);

F77_RET_I slaeda_f77(int *n, int *tlvls, int *curlvl, int *curpbm, int *prmptr, int *perm, int *givptr, int *givcol, float *givnum, float *q, int *qptr, float *z, float *ztemp, int *info);

F77_RET_I slaein_f77(logical *rightv, logical *noinit, int *n, float *h, int *ldh, float *wr, float *wi, float *vr, float *vi, float *b, int *ldb, float *work, float *eps3, float *smlnum, float *bignum, int *info);

F77_RET_I slaev2_f77(float *a, float *b, float *c, float *rt1, float *rt2, float *cs1, float *sn1);

F77_RET_I slaexc_f77(logical *wantq, int *n, float *t, int *ldt, float *q, int *ldq, int *j1, int *n1, int *n2, float *work, int *info);

F77_RET_I slag2_f77(float *a, int *lda, float *b, int *ldb, float *safmin, float *scale1, float *scale2, float *wr1, float *wr2, float *wi);

F77_RET_I slag2d_f77(int *m, int *n, float *sa, int *ldsa, double *a, int *lda, int *info);

F77_RET_I slags2_f77(logical *upper, float *a1, float *a2, float *a3, float *b1, float *b2, float *b3, float *csu, float *snu, float *csv, float *snv, float *csq, float *snq);

F77_RET_I slagtf_f77(int *n, float *a, float *lambda, float *b, float *c, float *tol, float *d, int *in, int *info);

F77_RET_I slagtm_f77(const char *trans, int *n, int *nrhs, float *alpha, float *dl, float *d, float *du, float *x, int *ldx, float *beta, float *b, int *ldb);

F77_RET_I slagts_f77(int *job, int *n, float *a, float *b, float *c, float *d, int *in, float *y, float *tol, int *info);

F77_RET_I slagv2_f77(float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *csl, float *snl, float *csr, float *snr);

F77_RET_I slahqr_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, float *h, int *ldh, float *wr, float *wi, int *iloz, int *ihiz, float *z, int *ldz, int *info);

F77_RET_I slahr2_f77(int *n, int *k, int *nb, float *a, int *lda, float *tau, float *t, int *ldt, float *y, int *ldy);

F77_RET_I slahrd_f77(int *n, int *k, int *nb, float *a, int *lda, float *tau, float *t, int *ldt, float *y, int *ldy);

F77_RET_I slaic1_f77(int *job, int *j, float *x, float *sest, float *w, float *gamma, float *sestpr, float *s, float *c);

logical slaisnan_f77(float *sin1, float *sin2);

F77_RET_I slaln2_f77(logical *ltrans, int *na, int *nw, float *smin, float *ca, float *a, int *lda, float *d1, float *d2, float *b, int *ldb, float *wr, float *wi, float *x, int *ldx, float *scale, float *xnorm, int *info);

F77_RET_I slals0_f77(int *icompq, int *nl, int *nr, int *sqre, int *nrhs, float *b, int *ldb, float *bx, int *ldbx, int *perm, int *givptr, int *givcol, int *ldgcol, float *givnum, int *ldgnum, float *poles, float *difl, float *difr, float *z, int *k, float *c, float *s, float *work, int *info);

F77_RET_I slalsa_f77(int *icompq, int *smlsiz, int *n, int *nrhs, float *b, int *ldb, float *bx, int *ldbx, float *u, int *ldu, float *vt, int *k, float *difl, float *difr, float *z, float *poles, int *givptr, int *givcol, int *ldgcol, int *perm, float *givnum, float *c, float *s, float *work, int *iwork, int *info);

F77_RET_I slalsd_f77(const char *uplo, int *smlsiz, int *n, int *nrhs, float *d, float *e, float *b, int *ldb, float *rcond, int *rank, float *work, int *iwork, int *info);

F77_RET_I slamrg_f77(int *n1, int *n2, float *a, int *strd1, int *strd2, int *index);

F77_RET_I slaneg_f77(int *n, float *d, float *lld, float *sigma, float *pivmin, int *r);

F77_RET_F slangb_f77(const char *norm, int *n, int *kl, int *ku, float *ab, int *ldab, float *work);

F77_RET_F slange_f77(const char *norm, int *m, int *n, float *a, int *lda, float *work);

F77_RET_F slangt_f77(const char *norm, int *n, float *dl, float *d, float *du);

F77_RET_F slanhs_f77(const char *norm, int *n, float *a, int *lda, float *work);

F77_RET_F slansb_f77(const char *norm, const char *uplo, int *n, int *k, float *ab, int *ldab, float *work);

F77_RET_F slansp_f77(const char *norm, const char *uplo, int *n, float *ap, float *work);

F77_RET_F slanst_f77(const char *norm, int *n, float *d, float *e);

F77_RET_F slansy_f77(const char *norm, const char *uplo, int *n, float *a, int *lda, float *work);

F77_RET_F slantb_f77(const char *norm, const char *uplo, const char *diag, int *n, int *k, float *ab, int *ldab, float *work);

F77_RET_F slantp_f77(const char *norm, const char *uplo, const char *diag, int *n, float *ap, float *work);

F77_RET_F slantr_f77(const char *norm, const char *uplo, const char *diag, int *m, int *n, float *a, int *lda, float *work);

F77_RET_I slanv2_f77(float *a, float *b, float *c, float *d, float *rt1r, float *rt1i, float *rt2r, float *rt2i, float *cs, float *sn);

F77_RET_I slapll_f77(int *n, float *x, int *incx, float *y, int *incy, float *ssmin);

F77_RET_I slapmt_f77(logical *forwrd, int *m, int *n, float *x, int *ldx, int *k);

F77_RET_F slapy2_f77(float *x, float *y);

F77_RET_F slapy3_f77(float *x, float *y, float *z);

F77_RET_I slaqgb_f77(int *m, int *n, int *kl, int *ku, float *ab, int *ldab, float *r, float *c, float *rowcnd, float *colcnd, float *amax, const char *equed);

F77_RET_I slaqge_f77(int *m, int *n, float *a, int *lda, float *r, float *c, float *rowcnd, float *colcnd, float *amax, const char *equed);

F77_RET_I slaqp2_f77(int *m, int *n, int *offset, float *a, int *lda, int *jpvt, float *tau, float *vn1, float *vn2, float *work);

F77_RET_I slaqps_f77(int *m, int *n, int *offset, int *nb, int *kb, float *a, int *lda, int *jpvt, float *tau, float *vn1, float *vn2, float *auxv, float *f, int *ldf);

F77_RET_I slaqr0_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, float *h, int *ldh, float *wr, float *wi, int *iloz, int *ihiz, float *z, int *ldz, float *work, int *lwork, int *info);

F77_RET_I slaqr1_f77(int *n, float *h, int *ldh, float *sr1, float *si1, float *sr2, float *si2, float *v);

F77_RET_I slaqr2_f77(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, float *h, int *ldh, int *iloz, int *ihiz, float *z, int *ldz, int *ns, int *nd, float *sr, float *si, float *v, int *ldv, int *nh, float *t, int *ldt, int *nv, float *wv, int *ldwv, float *work, int *lwork);

F77_RET_I slaqr3_f77(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, float *h, int *ldh, int *iloz, int *ihiz, float *z, int *ldz, int *ns, int *nd, float *sr, float *si, float *v, int *ldv, int *nh, float *t, int *ldt, int *nv, float *wv, int *ldwv, float *work, int *lwork);

F77_RET_I slaqr4_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, float *h, int *ldh, float *wr, float *wi, int *iloz, int *ihiz, float *z, int *ldz, float *work, int *lwork, int *info);

F77_RET_I slaqr5_f77(logical *wantt, logical *wantz, int *kacc22, int *n, int *ktop, int *kbot, int *nshfts, float *sr, float *si, float *h, int *ldh, int *iloz, int *ihiz, float *z, int *ldz, float *v, int *ldv, float *u, int *ldu, int *nv, float *wv, int *ldwv, int *nh, float *wh, int *ldwh);

F77_RET_I slaqsb_f77(const char *uplo, int *n, int *kd, float *ab, int *ldab, float *s, float *scond, float *amax, const char *equed);

F77_RET_I slaqsp_f77(const char *uplo, int *n, float *ap, float *s, float *scond, float *amax, const char *equed);

F77_RET_I slaqsy_f77(const char *uplo, int *n, float *a, int *lda, float *s, float *scond, float *amax, const char *equed);

F77_RET_I slaqtr_f77(logical *ltran, logical *lfloat, int *n, float *t, int *ldt, float *b, float *w, float *scale, float *x, float *work, int *info);

F77_RET_I slar1v_f77(int *n, int *b1, int *bn, float *lambda, float *d, float *l, float *ld, float *lld, float *pivmin, float *gaptol, float *z, logical *wantnc, int *negcnt, float *ztz, float *mingma, int *r, int *isuppz, float *nrminv, float *resid, float *rqcorr, float *work);

F77_RET_I slar2v_f77(int *n, float *x, float *y, float *z, int *incx, float *c, float *s, int *incc);

F77_RET_I slarf_f77(const char *side, int *m, int *n, float *v, int *incv, float *tau, float *c, int *ldc, float *work);

F77_RET_I slarfb_f77(const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n, int *k, float *v, int *ldv, float *t, int *ldt, float *c, int *ldc, float *work, int *ldwork);

F77_RET_I slarfg_f77(int *n, float *alpha, float *x, int *incx, float *tau);

F77_RET_I slarft_f77(const char *direct, const char *storev, int *n, int *k, float *v, int *ldv, float *tau, float *t, int *ldt);

F77_RET_I slarfx_f77(const char *side, int *m, int *n, float *v, float *tau, float *c, int *ldc, float *work);

F77_RET_I slargv_f77(int *n, float *x, int *incx, float *y, int *incy, float *c, int *incc);

F77_RET_I slarnv_f77(int *idist, int *iseed, int *n, float *x);

F77_RET_I slarra_f77(int *n, float *d, float *e, float *e2, float *spltol, float *tnrm, int *nsplit, int *isplit, int *info);

F77_RET_I slarrb_f77(int *n, float *d, float *lld, int *ifirst, int *ilast, float *rtol1, float *rtol2, int *offset, float *w, float *wgap, float *werr, float *work, int *iwork, float *pivmin, float *spdiam, int *twist, int *info);

F77_RET_I slarrc_f77(const char *jobt, int *n, float *vl, float *vu, float *d, float *e, float *pivmin, int *eigcnt, int *lcnt, int *rcnt, int *info);

F77_RET_I slarrd_f77(const char *range, const char *order, int *n, float *vl, float *vu, int *il, int *iu, float *gers, float *reltol, float *d, float *e, float *e2, float *pivmin, int *nsplit, int *isplit, int *m, float *w, float *werr, float *wl, float *wu, int *iblock, int *indexw, float *work, int *iwork, int *info);

F77_RET_I slarre_f77(const char *range, int *n, float *vl, float *vu, int *il, int *iu, float *d, float *e, float *e2, float *rtol1, float *rtol2, float *spltol, int *nsplit, int *isplit, int *m, float *w, float *werr, float *wgap, int *iblock, int *indexw, float *gers, float *pivmin, float *work, int *iwork, int *info);

F77_RET_I slarrf_f77(int *n, float *d, float *l, float *ld, int *clstrt, int *clend, float *w, float *wgap, float *werr, float *spdiam, float *clgapl, float *clgapr, float *pivmin, float *sigma, float *dplus, float *lplus, float *work, int *info);

F77_RET_I slarrj_f77(int *n, float *d, float *e2, int *ifirst, int *ilast, float *rtol, int *offset, float *w, float *werr, float *work, int *iwork, float *pivmin, float *spdiam, int *info);

F77_RET_I slarrk_f77(int *n, int *iw, float *gl, float *gu, float *d, float *e2, float *pivmin, float *reltol, float *w, float *werr, int *info);

F77_RET_I slarrr_f77(int *n, float *d, float *e, int *info);

F77_RET_I slarrv_f77(int *n, float *vl, float *vu, float *d, float *l, float *pivmin, int *isplit, int *m, int *dol, int *dou, float *minrgp, float *rtol1, float *rtol2, float *w, float *werr, float *wgap, int *iblock, int *indexw, float *gers, float *z, int *ldz, int *isuppz, float *work, int *iwork, int *info);

F77_RET_I slartg_f77(float *f, float *g, float *cs, float *sn, float *r);

F77_RET_I slartv_f77(int *n, float *x, int *incx, float *y, int *incy, float *c, float *s, int *incc);

F77_RET_I slaruv_f77(int *iseed, int *n, float *x);

F77_RET_I slarz_f77(const char *side, int *m, int *n, int *l, float *v, int *incv, float *tau, float *c, int *ldc, float *work);

F77_RET_I slarzb_f77(const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n, int *k, int *l, float *v, int *ldv, float *t, int *ldt, float *c, int *ldc, float *work, int *ldwork);

F77_RET_I slarzt_f77(const char *direct, const char *storev, int *n, int *k, float *v, int *ldv, float *tau, float *t, int *ldt);

F77_RET_I slas2_f77(float *f, float *g, float *h, float *ssmin, float *ssmax);

F77_RET_I slascl_f77(const char *type, int *kl, int *ku, float *cfrom, float *cto, int *m, int *n, float *a, int *lda, int *info);

F77_RET_I slasd0_f77(int *n, int *sqre, float *d, float *e, float *u, int *ldu, float *vt, int *ldvt, int *smlsiz, int *iwork, float *work, int *info);

F77_RET_I slasd1_f77(int *nl, int *nr, int *sqre, float *d, float *alpha, float *beta, float *u, int *ldu, float *vt, int *ldvt, int *idxq, int *iwork, float *work, int *info);

F77_RET_I slasd2_f77(int *nl, int *nr, int *sqre, int *k, float *d, float *z, float *alpha, float *beta, float *u, int *ldu, float *vt, int *ldvt, float *dsigma, float *u2, int *ldu2, float *vt2, int *ldvt2, int *idxp, int *idx, int *idxc, int *idxq, int *coltyp, int *info);

F77_RET_I slasd3_f77(int *nl, int *nr, int *sqre, int *k, float *d, float *q, int *ldq, float *dsigma, float *u, int *ldu, float *u2, int *ldu2, float *vt, int *ldvt, float *vt2, int *ldvt2, int *idxc, int *ctot, float *z, int *info);

F77_RET_I slasd4_f77(int *n, int *i, float *d, float *z, float *delta, float *rho, float *sigma, float *work, int *info);

F77_RET_I slasd5_f77(int *i, float *d, float *z, float *delta, float *rho, float *dsigma, float *work);

F77_RET_I slasd6_f77(int *icompq, int *nl, int *nr, int *sqre, float *d, float *vf, float *vl, float *alpha, float *beta, int *idxq, int *perm, int *givptr, int *givcol, int *ldgcol, float *givnum, int *ldgnum, float *poles, float *difl, float *difr, float *z, int *k, float *c, float *s, float *work, int *iwork, int *info);

F77_RET_I slasd7_f77(int *icompq, int *nl, int *nr, int *sqre, int *k, float *d, float *z, float *zw, float *vf, float *vfw, float *vl, float *vlw, float *alpha, float *beta, float *dsigma, int *idx, int *idxp, int *idxq, int *perm, int *givptr, int *givcol, int *ldgcol, float *givnum, int *ldgnum, float *c, float *s, int *info);

F77_RET_I slasd8_f77(int *icompq, int *k, float *d, float *z, float *vf, float *vl, float *difl, float *difr, int *lddifr, float *dsigma, float *work, int *info);

F77_RET_I slasda_f77(int *icompq, int *smlsiz, int *n, int *sqre, float *d, float *e, float *u, int *ldu, float *vt, int *k, float *difl, float *difr, float *z, float *poles, int *givptr, int *givcol, int *ldgcol, int *perm, float *givnum, float *c, float *s, float *work, int *iwork, int *info);

F77_RET_I slasdq_f77(const char *uplo, int *sqre, int *n, int *ncvt, int *nru, int *ncc, float *d, float *e, float *vt, int *ldvt, float *u, int *ldu, float *c, int *ldc, float *work, int *info);

F77_RET_I slasdt_f77(int *n, int *lvl, int *nd, int *inode, int *ndiml, int *ndimr, int *msub);

F77_RET_I slaset_f77(const char *uplo, int *m, int *n, float *alpha, float *beta, float *a, int *lda);

F77_RET_I slasq1_f77(int *n, float *d, float *e, float *work, int *info);

F77_RET_I slasq2_f77(int *n, float *z, int *info);

F77_RET_I slasq3_f77(int *i0, int *n0, float *z, int *pp, float *dmin, float *sigma, float *desig, float *qmax, int *nfail, int *iter, int *ndiv, logical *ieee);

F77_RET_I slasq4_f77(int *i0, int *n0, float *z, int *pp, int *n0in, float *dmin, float *dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *tau, int *ttype);

F77_RET_I slasq5_f77(int *i0, int *n0, float *z, int *pp, float *tau, float *dmin, float *dmin1, float *dmin2, float *dn, float *dnm1, float *dnm2, logical *ieee);

F77_RET_I slasq6_f77(int *i0, int *n0, float *z, int *pp, float *dmin, float *dmin1, float *dmin2, float *dn, float *dnm1, float *dnm2);

F77_RET_I slasr_f77(const char *side, const char *pivot, const char *direct, int *m, int *n, float *c, float *s, float *a, int *lda);

F77_RET_I slasrt_f77(const char *id, int *n, float *d, int *info);

F77_RET_I slassq_f77(int *n, float *x, int *incx, float *scale, float *sumsq);

F77_RET_I slasv2_f77(float *f, float *g, float *h, float *ssmin, float *ssmax, float *snr, float *csr, float *snl, float *csl);

F77_RET_I slaswp_f77(int *n, float *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);

F77_RET_I slasy2_f77(logical *ltranl, logical *ltranr, int *isgn, int *n1, int *n2, float *tl, int *ldtl, float *tr, int *ldtr, float *b, int *ldb, float *scale, float *x, int *ldx, float *xnorm, int *info);

F77_RET_I slasyf_f77(const char *uplo, int *n, int *nb, int *kb, float *a, int *lda, int *ipiv, float *w, int *ldw, int *info);

F77_RET_I slatbs_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, int *kd, float *ab, int *ldab, float *x, float *scale, float *cnorm, int *info);

F77_RET_I slatdf_f77(int *ijob, int *n, float *z, int *ldz, float *rhs, float *rdsum, float *rdscal, int *ipiv, int *jpiv);

F77_RET_I slatps_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, float *ap, float *x, float *scale, float *cnorm, int *info);

F77_RET_I slatrd_f77(const char *uplo, int *n, int *nb, float *a, int *lda, float *e, float *tau, float *w, int *ldw);

F77_RET_I slatrs_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, float *a, int *lda, float *x, float *scale, float *cnorm, int *info);

F77_RET_I slatrz_f77(int *m, int *n, int *l, float *a, int *lda, float *tau, float *work);

F77_RET_I slatzm_f77(const char *side, int *m, int *n, float *v, int *incv, float *tau, float *c1, float *c2, int *ldc, float *work);

F77_RET_I slauu2_f77(const char *uplo, int *n, float *a, int *lda, int *info);

F77_RET_I slauum_f77(const char *uplo, int *n, float *a, int *lda, int *info);

F77_RET_I slazq3_f77(int *i0, int *n0, float *z, int *pp, float *dmin, float *sigma, float *desig, float *qmax, int *nfail, int *iter, int *ndiv, logical *ieee, int *ttype, float *dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *tau);

F77_RET_I slazq4_f77(int *i0, int *n0, float *z, int *pp, int *n0in, float *dmin, float *dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *tau, int *ttype, float *g);

F77_RET_I sopgtr_f77(const char *uplo, int *n, float *ap, float *tau, float *q, int *ldq, float *work, int *info);

F77_RET_I sopmtr_f77(const char *side, const char *uplo, const char *trans, int *m, int *n, float *ap, float *tau, float *c, int *ldc, float *work, int *info);

F77_RET_I sorg2l_f77(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *info);

F77_RET_I sorg2r_f77(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *info);

F77_RET_I sorgbr_f77(const char *vect, int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sorghr_f77(int *n, int *ilo, int *ihi, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sorgl2_f77(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *info);

F77_RET_I sorglq_f77(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sorgql_f77(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sorgqr_f77(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sorgr2_f77(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *info);

F77_RET_I sorgrq_f77(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sorgtr_f77(const char *uplo, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I sorm2l_f77(const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *info);

F77_RET_I sorm2r_f77(const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *info);

F77_RET_I sormbr_f77(const char *vect, const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *lwork, int *info);

F77_RET_I sormhr_f77(const char *side, const char *trans, int *m, int *n, int *ilo, int *ihi, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *lwork, int *info);

F77_RET_I sorml2_f77(const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *info);

F77_RET_I sormlq_f77(const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *lwork, int *info);

F77_RET_I sormql_f77(const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *lwork, int *info);

F77_RET_I sormqr_f77(const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *lwork, int *info);

F77_RET_I sormr2_f77(const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *info);

F77_RET_I sormr3_f77(const char *side, const char *trans, int *m, int *n, int *k, int *l, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *info);

F77_RET_I sormrq_f77(const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *lwork, int *info);

F77_RET_I sormrz_f77(const char *side, const char *trans, int *m, int *n, int *k, int *l, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *lwork, int *info);

F77_RET_I sormtr_f77(const char *side, const char *uplo, const char *trans, int *m, int *n, float *a, int *lda, float *tau, float *c, int *ldc, float *work, int *lwork, int *info);

F77_RET_I spbcon_f77(const char *uplo, int *n, int *kd, float *ab, int *ldab, float *anorm, float *rcond, float *work, int *iwork, int *info);

F77_RET_I spbequ_f77(const char *uplo, int *n, int *kd, float *ab, int *ldab, float *s, float *scond, float *amax, int *info);

F77_RET_I spbrfs_f77(const char *uplo, int *n, int *kd, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I spbstf_f77(const char *uplo, int *n, int *kd, float *ab, int *ldab, int *info);

F77_RET_I spbsv_f77(const char *uplo, int *n, int *kd, int *nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);

F77_RET_I spbsvx_f77(const char *fact, const char *uplo, int *n, int *kd, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, const char *equed, float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I spbtf2_f77(const char *uplo, int *n, int *kd, float *ab, int *ldab, int *info);

F77_RET_I spbtrf_f77(const char *uplo, int *n, int *kd, float *ab, int *ldab, int *info);

F77_RET_I spbtrs_f77(const char *uplo, int *n, int *kd, int *nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);

F77_RET_I spocon_f77(const char *uplo, int *n, float *a, int *lda, float *anorm, float *rcond, float *work, int *iwork, int *info);

F77_RET_I spoequ_f77(int *n, float *a, int *lda, float *s, float *scond, float *amax, int *info);

F77_RET_I sporfs_f77(const char *uplo, int *n, int *nrhs, float *a, int *lda, float *af, int *ldaf, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I sposv_f77(const char *uplo, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, int *info);

F77_RET_I sposvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, float *a, int *lda, float *af, int *ldaf, const char *equed, float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I spotf2_f77(const char *uplo, int *n, float *a, int *lda, int *info);

F77_RET_I spotrf_f77(const char *uplo, int *n, float *a, int *lda, int *info);

F77_RET_I spotri_f77(const char *uplo, int *n, float *a, int *lda, int *info);

F77_RET_I spotrs_f77(const char *uplo, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, int *info);

F77_RET_I sppcon_f77(const char *uplo, int *n, float *ap, float *anorm, float *rcond, float *work, int *iwork, int *info);

F77_RET_I sppequ_f77(const char *uplo, int *n, float *ap, float *s, float *scond, float *amax, int *info);

F77_RET_I spprfs_f77(const char *uplo, int *n, int *nrhs, float *ap, float *afp, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I sppsv_f77(const char *uplo, int *n, int *nrhs, float *ap, float *b, int *ldb, int *info);

F77_RET_I sppsvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, float *ap, float *afp, const char *equed, float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I spptrf_f77(const char *uplo, int *n, float *ap, int *info);

F77_RET_I spptri_f77(const char *uplo, int *n, float *ap, int *info);

F77_RET_I spptrs_f77(const char *uplo, int *n, int *nrhs, float *ap, float *b, int *ldb, int *info);

F77_RET_I sptcon_f77(int *n, float *d, float *e, float *anorm, float *rcond, float *work, int *info);

F77_RET_I spteqr_f77(const char *compz, int *n, float *d, float *e, float *z, int *ldz, float *work, int *info);

F77_RET_I sptrfs_f77(int *n, int *nrhs, float *d, float *e, float *df, float *ef, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *info);

F77_RET_I sptsv_f77(int *n, int *nrhs, float *d, float *e, float *b, int *ldb, int *info);

F77_RET_I sptsvx_f77(const char *fact, int *n, int *nrhs, float *d, float *e, float *df, float *ef, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *info);

F77_RET_I spttrf_f77(int *n, float *d, float *e, int *info);

F77_RET_I spttrs_f77(int *n, int *nrhs, float *d, float *e, float *b, int *ldb, int *info);

F77_RET_I sptts2_f77(int *n, int *nrhs, float *d, float *e, float *b, int *ldb);

F77_RET_I srscl_f77(int *n, float *sa, float *sx, int *incx);

F77_RET_I ssbev_f77(const char *jobz, const char *uplo, int *n, int *kd, float *ab, int *ldab, float *w, float *z, int *ldz, float *work, int *info);

F77_RET_I ssbevd_f77(const char *jobz, const char *uplo, int *n, int *kd, float *ab, int *ldab, float *w, float *z, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I ssbevx_f77(const char *jobz, const char *range, const char *uplo, int *n, int *kd, float *ab, int *ldab, float *q, int *ldq, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z, int *ldz, float *work, int *iwork, int *ifail, int *info);

F77_RET_I ssbgst_f77(const char *vect, const char *uplo, int *n, int *ka, int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *x, int *ldx, float *work, int *info);

F77_RET_I ssbgv_f77(const char *jobz, const char *uplo, int *n, int *ka, int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *w, float *z, int *ldz, float *work, int *info);

F77_RET_I ssbgvd_f77(const char *jobz, const char *uplo, int *n, int *ka, int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *w, float *z, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I ssbgvx_f77(const char *jobz, const char *range, const char *uplo, int *n, int *ka, int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *q, int *ldq, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z, int *ldz, float *work, int *iwork, int *ifail, int *info);

F77_RET_I ssbtrd_f77(const char *vect, const char *uplo, int *n, int *kd, float *ab, int *ldab, float *d, float *e, float *q, int *ldq, float *work, int *info);

F77_RET_I sspcon_f77(const char *uplo, int *n, float *ap, int *ipiv, float *anorm, float *rcond, float *work, int *iwork, int *info);

F77_RET_I sspev_f77(const char *jobz, const char *uplo, int *n, float *ap, float *w, float *z, int *ldz, float *work, int *info);

F77_RET_I sspevd_f77(const char *jobz, const char *uplo, int *n, float *ap, float *w, float *z, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I sspevx_f77(const char *jobz, const char *range, const char *uplo, int *n, float *ap, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z, int *ldz, float *work, int *iwork, int *ifail, int *info);

F77_RET_I sspgst_f77(int *itype, const char *uplo, int *n, float *ap, float *bp, int *info);

F77_RET_I sspgv_f77(int *itype, const char *jobz, const char *uplo, int *n, float *ap, float *bp, float *w, float *z, int *ldz, float *work, int *info);

F77_RET_I sspgvd_f77(int *itype, const char *jobz, const char *uplo, int *n, float *ap, float *bp, float *w, float *z, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I sspgvx_f77(int *itype, const char *jobz, const char *range, const char *uplo, int *n, float *ap, float *bp, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z, int *ldz, float *work, int *iwork, int *ifail, int *info);

F77_RET_I ssprfs_f77(const char *uplo, int *n, int *nrhs, float *ap, float *afp, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I sspsv_f77(const char *uplo, int *n, int *nrhs, float *ap, int *ipiv, float *b, int *ldb, int *info);

F77_RET_I sspsvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, float *ap, float *afp, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I ssptrd_f77(const char *uplo, int *n, float *ap, float *d, float *e, float *tau, int *info);

F77_RET_I ssptrf_f77(const char *uplo, int *n, float *ap, int *ipiv, int *info);

F77_RET_I ssptri_f77(const char *uplo, int *n, float *ap, int *ipiv, float *work, int *info);

F77_RET_I ssptrs_f77(const char *uplo, int *n, int *nrhs, float *ap, int *ipiv, float *b, int *ldb, int *info);

F77_RET_I sstebz_f77(const char *range, const char *order, int *n, float *vl, float *vu, int *il, int *iu, float *abstol, float *d, float *e, int *m, int *nsplit, float *w, int *iblock, int *isplit, float *work, int *iwork, int *info);

F77_RET_I sstedc_f77(const char *compz, int *n, float *d, float *e, float *z, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I sstegr_f77(const char *jobz, const char *range, int *n, float *d, float *e, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z, int *ldz, int *isuppz, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I sstein_f77(int *n, float *d, float *e, int *m, float *w, int *iblock, int *isplit, float *z, int *ldz, float *work, int *iwork, int *ifail, int *info);

F77_RET_I sstemr_f77(const char *jobz, const char *range, int *n, float *d, float *e, float *vl, float *vu, int *il, int *iu, int *m, float *w, float *z, int *ldz, int *nzc, int *isuppz, logical *tryrac, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I ssteqr_f77(const char *compz, int *n, float *d, float *e, float *z, int *ldz, float *work, int *info);

F77_RET_I ssterf_f77(int *n, float *d, float *e, int *info);

F77_RET_I sstev_f77(const char *jobz, int *n, float *d, float *e, float *z, int *ldz, float *work, int *info);

F77_RET_I sstevd_f77(const char *jobz, int *n, float *d, float *e, float *z, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I sstevr_f77(const char *jobz, const char *range, int *n, float *d, float *e, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z, int *ldz, int *isuppz, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I sstevx_f77(const char *jobz, const char *range, int *n, float *d, float *e, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z, int *ldz, float *work, int *iwork, int *ifail, int *info);

F77_RET_I ssycon_f77(const char *uplo, int *n, float *a, int *lda, int *ipiv, float *anorm, float *rcond, float *work, int *iwork, int *info);

F77_RET_I ssyev_f77(const char *jobz, const char *uplo, int *n, float *a, int *lda, float *w, float *work, int *lwork, int *info);

F77_RET_I ssyevd_f77(const char *jobz, const char *uplo, int *n, float *a, int *lda, float *w, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I ssyevr_f77(const char *jobz, const char *range, const char *uplo, int *n, float *a, int *lda, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z, int *ldz, int *isuppz, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I ssyevx_f77(const char *jobz, const char *range, const char *uplo, int *n, float *a, int *lda, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z, int *ldz, float *work, int *lwork, int *iwork, int *ifail, int *info);

F77_RET_I ssygs2_f77(int *itype, const char *uplo, int *n, float *a, int *lda, float *b, int *ldb, int *info);

F77_RET_I ssygst_f77(int *itype, const char *uplo, int *n, float *a, int *lda, float *b, int *ldb, int *info);

F77_RET_I ssygv_f77(int *itype, const char *jobz, const char *uplo, int *n, float *a, int *lda, float *b, int *ldb, float *w, float *work, int *lwork, int *info);

F77_RET_I ssygvd_f77(int *itype, const char *jobz, const char *uplo, int *n, float *a, int *lda, float *b, int *ldb, float *w, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I ssygvx_f77(int *itype, const char *jobz, const char *range, const char *uplo, int *n, float *a, int *lda, float *b, int *ldb, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z, int *ldz, float *work, int *lwork, int *iwork, int *ifail, int *info);

F77_RET_I ssyrfs_f77(const char *uplo, int *n, int *nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I ssysv_f77(const char *uplo, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, float *work, int *lwork, int *info);

F77_RET_I ssysvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *lwork, int *iwork, int *info);

F77_RET_I ssytd2_f77(const char *uplo, int *n, float *a, int *lda, float *d, float *e, float *tau, int *info);

F77_RET_I ssytf2_f77(const char *uplo, int *n, float *a, int *lda, int *ipiv, int *info);

F77_RET_I ssytrd_f77(const char *uplo, int *n, float *a, int *lda, float *d, float *e, float *tau, float *work, int *lwork, int *info);

F77_RET_I ssytrf_f77(const char *uplo, int *n, float *a, int *lda, int *ipiv, float *work, int *lwork, int *info);

F77_RET_I ssytri_f77(const char *uplo, int *n, float *a, int *lda, int *ipiv, float *work, int *info);

F77_RET_I ssytrs_f77(const char *uplo, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);

F77_RET_I stbcon_f77(const char *norm, const char *uplo, const char *diag, int *n, int *kd, float *ab, int *ldab, float *rcond, float *work, int *iwork, int *info);

F77_RET_I stbrfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *kd, int *nrhs, float *ab, int *ldab, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I stbtrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *kd, int *nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);

F77_RET_I stgevc_f77(const char *side, const char *howmny, logical *select, int *n, float *s, int *lds, float *p, int *ldp, float *vl, int *ldvl, float *vr, int *ldvr, int *mm, int *m, float *work, int *info);

F77_RET_I stgex2_f77(logical *wantq, logical *wantz, int *n, float *a, int *lda, float *b, int *ldb, float *q, int *ldq, float *z, int *ldz, int *j1, int *n1, int *n2, float *work, int *lwork, int *info);

F77_RET_I stgexc_f77(logical *wantq, logical *wantz, int *n, float *a, int *lda, float *b, int *ldb, float *q, int *ldq, float *z, int *ldz, int *ifst, int *ilst, float *work, int *lwork, int *info);

F77_RET_I stgsen_f77(int *ijob, logical *wantq, logical *wantz, logical *select, int *n, float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *q, int *ldq, float *z, int *ldz, int *m, float *pl, float *pr, float *dif, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I stgsja_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *p, int *n, int *k, int *l, float *a, int *lda, float *b, int *ldb, float *tola, float *tolb, float *alpha, float *beta, float *u, int *ldu, float *v, int *ldv, float *q, int *ldq, float *work, int *ncycle, int *info);

F77_RET_I stgsna_f77(const char *job, const char *howmny, logical *select, int *n, float *a, int *lda, float *b, int *ldb, float *vl, int *ldvl, float *vr, int *ldvr, float *s, float *dif, int *mm, int *m, float *work, int *lwork, int *iwork, int *info);

F77_RET_I stgsy2_f77(const char *trans, int *ijob, int *m, int *n, float *a, int *lda, float *b, int *ldb, float *c, int *ldc, float *d, int *ldd, float *e, int *lde, float *f, int *ldf, float *scale, float *rdsum, float *rdscal, int *iwork, int *pq, int *info);

F77_RET_I stgsyl_f77(const char *trans, int *ijob, int *m, int *n, float *a, int *lda, float *b, int *ldb, float *c, int *ldc, float *d, int *ldd, float *e, int *lde, float *f, int *ldf, float *scale, float *dif, float *work, int *lwork, int *iwork, int *info);

F77_RET_I stpcon_f77(const char *norm, const char *uplo, const char *diag, int *n, float *ap, float *rcond, float *work, int *iwork, int *info);

F77_RET_I stprfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, float *ap, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I stptri_f77(const char *uplo, const char *diag, int *n, float *ap, int *info);

F77_RET_I stptrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, float *ap, float *b, int *ldb, int *info);

F77_RET_I strcon_f77(const char *norm, const char *uplo, const char *diag, int *n, float *a, int *lda, float *rcond, float *work, int *iwork, int *info);

F77_RET_I strevc_f77(const char *side, const char *howmny, logical *select, int *n, float *t, int *ldt, float *vl, int *ldvl, float *vr, int *ldvr, int *mm, int *m, float *work, int *info);

F77_RET_I strexc_f77(const char *compq, int *n, float *t, int *ldt, float *q, int *ldq, int *ifst, int *ilst, float *work, int *info);

F77_RET_I strrfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);

F77_RET_I strsen_f77(const char *job, const char *compq, logical *select, int *n, float *t, int *ldt, float *q, int *ldq, float *wr, float *wi, int *m, float *s, float *sep, float *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I strsna_f77(const char *job, const char *howmny, logical *select, int *n, float *t, int *ldt, float *vl, int *ldvl, float *vr, int *ldvr, float *s, float *sep, int *mm, int *m, float *work, int *ldwork, int *iwork, int *info);

F77_RET_I strsyl_f77(const char *trana, const char *tranb, int *isgn, int *m, int *n, float *a, int *lda, float *b, int *ldb, float *c, int *ldc, float *scale, int *info);

F77_RET_I strti2_f77(const char *uplo, const char *diag, int *n, float *a, int *lda, int *info);

F77_RET_I strtri_f77(const char *uplo, const char *diag, int *n, float *a, int *lda, int *info);

F77_RET_I strtrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, int *info);

F77_RET_I stzrqf_f77(int *m, int *n, float *a, int *lda, float *tau, int *info);

F77_RET_I stzrzf_f77(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

F77_RET_I xerbla_f77(const char *srname, int *info);

F77_RET_I zbdsqr_f77(const char *uplo, int *n, int *ncvt, int *nru, int *ncc, double *d, double *e, doublecomplex *vt, int *ldvt, doublecomplex *u, int *ldu, doublecomplex *c, int *ldc, double *rwork, int *info);

F77_RET_I zcgesv_f77(int *n, int *nrhs, doublecomplex *a, int *lda, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, doublecomplex *work, singlecomplex *swork, int *iter, int *info);

F77_RET_I zdrscl_f77(int *n, double *sa, doublecomplex *sx, int *incx);

F77_RET_I zgbbrd_f77(const char *vect, int *m, int *n, int *ncc, int *kl, int *ku, doublecomplex *ab, int *ldab, double *d, double *e, doublecomplex *q, int *ldq, doublecomplex *pt, int *ldpt, doublecomplex *c, int *ldc, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgbcon_f77(const char *norm, int *n, int *kl, int *ku, doublecomplex *ab, int *ldab, int *ipiv, double *anorm, double *rcond, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgbequ_f77(int *m, int *n, int *kl, int *ku, doublecomplex *ab, int *ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info);

F77_RET_I zgbrfs_f77(const char *trans, int *n, int *kl, int *ku, int *nrhs, doublecomplex *ab, int *ldab, doublecomplex *afb, int *ldafb, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgbsv_f77(int *n, int *kl, int *ku, int *nrhs, doublecomplex *ab, int *ldab, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I zgbsvx_f77(const char *fact, const char *trans, int *n, int *kl, int *ku, int *nrhs, doublecomplex *ab, int *ldab, doublecomplex *afb, int *ldafb, int *ipiv, const char *equed, double *r, double *c, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgbtf2_f77(int *m, int *n, int *kl, int *ku, doublecomplex *ab, int *ldab, int *ipiv, int *info);

F77_RET_I zgbtrf_f77(int *m, int *n, int *kl, int *ku, doublecomplex *ab, int *ldab, int *ipiv, int *info);

F77_RET_I zgbtrs_f77(const char *trans, int *n, int *kl, int *ku, int *nrhs, doublecomplex *ab, int *ldab, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I zgebak_f77(const char *job, const char *side, int *n, int *ilo, int *ihi, double *scale, int *m, doublecomplex *v, int *ldv, int *info);

F77_RET_I zgebal_f77(const char *job, int *n, doublecomplex *a, int *lda, int *ilo, int *ihi, double *scale, int *info);

F77_RET_I zgebd2_f77(int *m, int *n, doublecomplex *a, int *lda, double *d, double *e, doublecomplex *tauq, doublecomplex *taup, doublecomplex *work, int *info);

F77_RET_I zgebrd_f77(int *m, int *n, doublecomplex *a, int *lda, double *d, double *e, doublecomplex *tauq, doublecomplex *taup, doublecomplex *work, int *lwork, int *info);

F77_RET_I zgecon_f77(const char *norm, int *n, doublecomplex *a, int *lda, double *anorm, double *rcond, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgeequ_f77(int *m, int *n, doublecomplex *a, int *lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info);

F77_RET_I zgees_f77(const char *jobvs, const char *sort, L_fp select, int *n, doublecomplex *a, int *lda, int *sdim, doublecomplex *w, doublecomplex *vs, int *ldvs, doublecomplex *work, int *lwork, double *rwork, logical *bwork, int *info);

F77_RET_I zgeesx_f77(const char *jobvs, const char *sort, L_fp select, const char *sense, int *n, doublecomplex *a, int *lda, int *sdim, doublecomplex *w, doublecomplex *vs, int *ldvs, double *rconde, double *rcondv, doublecomplex *work, int *lwork, double *rwork, logical *bwork, int *info);

F77_RET_I zgeev_f77(const char *jobvl, const char *jobvr, int *n, doublecomplex *a, int *lda, doublecomplex *w, doublecomplex *vl, int *ldvl, doublecomplex *vr, int *ldvr, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zgeevx_f77(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, int *n, doublecomplex *a, int *lda, doublecomplex *w, doublecomplex *vl, int *ldvl, doublecomplex *vr, int *ldvr, int *ilo, int *ihi, double *scale, double *abnrm, double *rconde, double *rcondv, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zgegs_f77(const char *jobvsl, const char *jobvsr, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, int *ldvsl, doublecomplex *vsr, int *ldvsr, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zgegv_f77(const char *jobvl, const char *jobvr, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, int *ldvl, doublecomplex *vr, int *ldvr, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zgehd2_f77(int *n, int *ilo, int *ihi, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *info);

F77_RET_I zgehrd_f77(int *n, int *ilo, int *ihi, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zgelq2_f77(int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *info);

F77_RET_I zgelqf_f77(int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zgels_f77(const char *trans, int *m, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *work, int *lwork, int *info);

F77_RET_I zgelsd_f77(int *m, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, double *s, double *rcond, int *rank, doublecomplex *work, int *lwork, double *rwork, int *iwork, int *info);

F77_RET_I zgelss_f77(int *m, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, double *s, double *rcond, int *rank, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zgelsx_f77(int *m, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, int *jpvt, double *rcond, int *rank, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgelsy_f77(int *m, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, int *jpvt, double *rcond, int *rank, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zgeql2_f77(int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *info);

F77_RET_I zgeqlf_f77(int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zgeqp3_f77(int *m, int *n, doublecomplex *a, int *lda, int *jpvt, doublecomplex *tau, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zgeqpf_f77(int *m, int *n, doublecomplex *a, int *lda, int *jpvt, doublecomplex *tau, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgeqr2_f77(int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *info);

F77_RET_I zgeqrf_f77(int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zgerfs_f77(const char *trans, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *af, int *ldaf, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgerq2_f77(int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *info);

F77_RET_I zgerqf_f77(int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zgesc2_f77(int *n, doublecomplex *a, int *lda, doublecomplex *rhs, int *ipiv, int *jpiv, double *scale);

F77_RET_I zgesdd_f77(const char *jobz, int *m, int *n, doublecomplex *a, int *lda, double *s, doublecomplex *u, int *ldu, doublecomplex *vt, int *ldvt, doublecomplex *work, int *lwork, double *rwork, int *iwork, int *info);

F77_RET_I zgesv_f77(int *n, int *nrhs, doublecomplex *a, int *lda, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I zgesvd_f77(const char *jobu, const char *jobvt, int *m, int *n, doublecomplex *a, int *lda, double *s, doublecomplex *u, int *ldu, doublecomplex *vt, int *ldvt, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zgesvx_f77(const char *fact, const char *trans, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *af, int *ldaf, int *ipiv, const char *equed, double *r, double *c, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgetc2_f77(int *n, doublecomplex *a, int *lda, int *ipiv, int *jpiv, int *info);

F77_RET_I zgetf2_f77(int *m, int *n, doublecomplex *a, int *lda, int *ipiv, int *info);

F77_RET_I zgetrf_f77(int *m, int *n, doublecomplex *a, int *lda, int *ipiv, int *info);

F77_RET_I zgetri_f77(int *n, doublecomplex *a, int *lda, int *ipiv, doublecomplex *work, int *lwork, int *info);

F77_RET_I zgetrs_f77(const char *trans, int *n, int *nrhs, doublecomplex *a, int *lda, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I zggbak_f77(const char *job, const char *side, int *n, int *ilo, int *ihi, double *lscale, double *rscale, int *m, doublecomplex *v, int *ldv, int *info);

F77_RET_I zggbal_f77(const char *job, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, int *ilo, int *ihi, double *lscale, double *rscale, double *work, int *info);

F77_RET_I zgges_f77(const char *jobvsl, const char *jobvsr, const char *sort, L_fp selctg, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, int *sdim, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, int *ldvsl, doublecomplex *vsr, int *ldvsr, doublecomplex *work, int *lwork, double *rwork, logical *bwork, int *info);

F77_RET_I zggesx_f77(const char *jobvsl, const char *jobvsr, const char *sort, L_fp selctg, const char *sense, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, int *sdim, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, int *ldvsl, doublecomplex *vsr, int *ldvsr, double *rconde, double *rcondv, doublecomplex *work, int *lwork, double *rwork, int *iwork, int *liwork, logical *bwork, int *info);

F77_RET_I zggev_f77(const char *jobvl, const char *jobvr, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, int *ldvl, doublecomplex *vr, int *ldvr, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zggevx_f77(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, int *ldvl, doublecomplex *vr, int *ldvr, int *ilo, int *ihi, double *lscale, double *rscale, double *abnrm, double *bbnrm, double *rconde, double *rcondv, doublecomplex *work, int *lwork, double *rwork, int *iwork, logical *bwork, int *info);

F77_RET_I zggglm_f77(int *n, int *m, int *p, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *d, doublecomplex *x, doublecomplex *y, doublecomplex *work, int *lwork, int *info);

F77_RET_I zgghrd_f77(const char *compq, const char *compz, int *n, int *ilo, int *ihi, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *q, int *ldq, doublecomplex *z, int *ldz, int *info);

F77_RET_I zgglse_f77(int *m, int *n, int *p, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *c, doublecomplex *d, doublecomplex *x, doublecomplex *work, int *lwork, int *info);

F77_RET_I zggqrf_f77(int *n, int *m, int *p, doublecomplex *a, int *lda, doublecomplex *taua, doublecomplex *b, int *ldb, doublecomplex *taub, doublecomplex *work, int *lwork, int *info);

F77_RET_I zggrqf_f77(int *m, int *p, int *n, doublecomplex *a, int *lda, doublecomplex *taua, doublecomplex *b, int *ldb, doublecomplex *taub, doublecomplex *work, int *lwork, int *info);

F77_RET_I zggsvd_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *n, int *p, int *k, int *l, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, double *alpha, double *beta, doublecomplex *u, int *ldu, doublecomplex *v, int *ldv, doublecomplex *q, int *ldq, doublecomplex *work, double *rwork, int *iwork, int *info);

F77_RET_I zggsvp_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *p, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, double *tola, double *tolb, int *k, int *l, doublecomplex *u, int *ldu, doublecomplex *v, int *ldv, doublecomplex *q, int *ldq, int *iwork, double *rwork, doublecomplex *tau, doublecomplex *work, int *info);

F77_RET_I zgtcon_f77(const char *norm, int *n, doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex *du2, int *ipiv, double *anorm, double *rcond, doublecomplex *work, int *info);

F77_RET_I zgtrfs_f77(const char *trans, int *n, int *nrhs, doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, doublecomplex *du2, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgtsv_f77(int *n, int *nrhs, doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex *b, int *ldb, int *info);

F77_RET_I zgtsvx_f77(const char *fact, const char *trans, int *n, int *nrhs, doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, doublecomplex *du2, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zgttrf_f77(int *n, doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex *du2, int *ipiv, int *info);

F77_RET_I zgttrs_f77(const char *trans, int *n, int *nrhs, doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex *du2, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I zgtts2_f77(int *itrans, int *n, int *nrhs, doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex *du2, int *ipiv, doublecomplex *b, int *ldb);

F77_RET_I zhbev_f77(const char *jobz, const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, double *w, doublecomplex *z, int *ldz, doublecomplex *work, double *rwork, int *info);

F77_RET_I zhbevd_f77(const char *jobz, const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, double *w, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I zhbevx_f77(const char *jobz, const char *range, const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, doublecomplex *q, int *ldq, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, doublecomplex *z, int *ldz, doublecomplex *work, double *rwork, int *iwork, int *ifail, int *info);

F77_RET_I zhbgst_f77(const char *vect, const char *uplo, int *n, int *ka, int *kb, doublecomplex *ab, int *ldab, doublecomplex *bb, int *ldbb, doublecomplex *x, int *ldx, doublecomplex *work, double *rwork, int *info);

F77_RET_I zhbgv_f77(const char *jobz, const char *uplo, int *n, int *ka, int *kb, doublecomplex *ab, int *ldab, doublecomplex *bb, int *ldbb, double *w, doublecomplex *z, int *ldz, doublecomplex *work, double *rwork, int *info);

F77_RET_I zhbgvd_f77(const char *jobz, const char *uplo, int *n, int *ka, int *kb, doublecomplex *ab, int *ldab, doublecomplex *bb, int *ldbb, double *w, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I zhbgvx_f77(const char *jobz, const char *range, const char *uplo, int *n, int *ka, int *kb, doublecomplex *ab, int *ldab, doublecomplex *bb, int *ldbb, doublecomplex *q, int *ldq, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, doublecomplex *z, int *ldz, doublecomplex *work, double *rwork, int *iwork, int *ifail, int *info);

F77_RET_I zhbtrd_f77(const char *vect, const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, double *d, double *e, doublecomplex *q, int *ldq, doublecomplex *work, int *info);

F77_RET_I zhecon_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *ipiv, double *anorm, double *rcond, doublecomplex *work, int *info);

F77_RET_I zheev_f77(const char *jobz, const char *uplo, int *n, doublecomplex *a, int *lda, double *w, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zheevd_f77(const char *jobz, const char *uplo, int *n, doublecomplex *a, int *lda, double *w, doublecomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I zheevr_f77(const char *jobz, const char *range, const char *uplo, int *n, doublecomplex *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, doublecomplex *z, int *ldz, int *isuppz, doublecomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I zheevx_f77(const char *jobz, const char *range, const char *uplo, int *n, doublecomplex *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, double *rwork, int *iwork, int *ifail, int *info);

F77_RET_I zhegs2_f77(int *itype, const char *uplo, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, int *info);

F77_RET_I zhegst_f77(int *itype, const char *uplo, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, int *info);

F77_RET_I zhegv_f77(int *itype, const char *jobz, const char *uplo, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, double *w, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zhegvd_f77(int *itype, const char *jobz, const char *uplo, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, double *w, doublecomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I zhegvx_f77(int *itype, const char *jobz, const char *range, const char *uplo, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, double *rwork, int *iwork, int *ifail, int *info);

F77_RET_I zherfs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *af, int *ldaf, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zhesv_f77(const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *work, int *lwork, int *info);

F77_RET_I zhesvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *af, int *ldaf, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zhetd2_f77(const char *uplo, int *n, doublecomplex *a, int *lda, double *d, double *e, doublecomplex *tau, int *info);

F77_RET_I zhetf2_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *ipiv, int *info);

F77_RET_I zhetrd_f77(const char *uplo, int *n, doublecomplex *a, int *lda, double *d, double *e, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zhetrf_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *ipiv, doublecomplex *work, int *lwork, int *info);

F77_RET_I zhetri_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *ipiv, doublecomplex *work, int *info);

F77_RET_I zhetrs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I zhgeqz_f77(const char *job, const char *compq, const char *compz, int *n, int *ilo, int *ihi, doublecomplex *h, int *ldh, doublecomplex *t, int *ldt, doublecomplex *alpha, doublecomplex *beta, doublecomplex *q, int *ldq, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zhpcon_f77(const char *uplo, int *n, doublecomplex *ap, int *ipiv, double *anorm, double *rcond, doublecomplex *work, int *info);

F77_RET_I zhpev_f77(const char *jobz, const char *uplo, int *n, doublecomplex *ap, double *w, doublecomplex *z, int *ldz, doublecomplex *work, double *rwork, int *info);

F77_RET_I zhpevd_f77(const char *jobz, const char *uplo, int *n, doublecomplex *ap, double *w, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I zhpevx_f77(const char *jobz, const char *range, const char *uplo, int *n, doublecomplex *ap, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, doublecomplex *z, int *ldz, doublecomplex *work, double *rwork, int *iwork, int *ifail, int *info);

F77_RET_I zhpgst_f77(int *itype, const char *uplo, int *n, doublecomplex *ap, doublecomplex *bp, int *info);

F77_RET_I zhpgv_f77(int *itype, const char *jobz, const char *uplo, int *n, doublecomplex *ap, doublecomplex *bp, double *w, doublecomplex *z, int *ldz, doublecomplex *work, double *rwork, int *info);

F77_RET_I zhpgvd_f77(int *itype, const char *jobz, const char *uplo, int *n, doublecomplex *ap, doublecomplex *bp, double *w, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I zhpgvx_f77(int *itype, const char *jobz, const char *range, const char *uplo, int *n, doublecomplex *ap, doublecomplex *bp, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, doublecomplex *z, int *ldz, doublecomplex *work, double *rwork, int *iwork, int *ifail, int *info);

F77_RET_I zhprfs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *ap, doublecomplex *afp, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zhpsv_f77(const char *uplo, int *n, int *nrhs, doublecomplex *ap, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I zhpsvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, doublecomplex *ap, doublecomplex *afp, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zhptrd_f77(const char *uplo, int *n, doublecomplex *ap, double *d, double *e, doublecomplex *tau, int *info);

F77_RET_I zhptrf_f77(const char *uplo, int *n, doublecomplex *ap, int *ipiv, int *info);

F77_RET_I zhptri_f77(const char *uplo, int *n, doublecomplex *ap, int *ipiv, doublecomplex *work, int *info);

F77_RET_I zhptrs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *ap, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I zhsein_f77(const char *side, const char *eigsrc, const char *initv, logical *select, int *n, doublecomplex *h, int *ldh, doublecomplex *w, doublecomplex *vl, int *ldvl, doublecomplex *vr, int *ldvr, int *mm, int *m, doublecomplex *work, double *rwork, int *ifaill, int *ifailr, int *info);

F77_RET_I zhseqr_f77(const char *job, const char *compz, int *n, int *ilo, int *ihi, doublecomplex *h, int *ldh, doublecomplex *w, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, int *info);

F77_RET_I zlabrd_f77(int *m, int *n, int *nb, doublecomplex *a, int *lda, double *d, double *e, doublecomplex *tauq, doublecomplex *taup, doublecomplex *x, int *ldx, doublecomplex *y, int *ldy);

F77_RET_I zlacgv_f77(int *n, doublecomplex *x, int *incx);

F77_RET_I zlacn2_f77(int *n, doublecomplex *v, doublecomplex *x, double *est, int *kase, int *isave);

F77_RET_I zlacon_f77(int *n, doublecomplex *v, doublecomplex *x, double *est, int *kase);

F77_RET_I zlacp2_f77(const char *uplo, int *m, int *n, double *a, int *lda, doublecomplex *b, int *ldb);

F77_RET_I zlacpy_f77(const char *uplo, int *m, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb);

F77_RET_I zlacrm_f77(int *m, int *n, doublecomplex *a, int *lda, double *b, int *ldb, doublecomplex *c, int *ldc, double *rwork);

F77_RET_I zlacrt_f77(int *n, doublecomplex *cx, int *incx, doublecomplex *cy, int *incy, doublecomplex *c, doublecomplex *s);

F77_RET_Z zladiv_f77(doublecomplex *x, doublecomplex *y);

F77_RET_I zlaed0_f77(int *qsiz, int *n, double *d, double *e, doublecomplex *q, int *ldq, doublecomplex *qstore, int *ldqs, double *rwork, int *iwork, int *info);

F77_RET_I zlaed7_f77(int *n, int *cutpnt, int *qsiz, int *tlvls, int *curlvl, int *curpbm, double *d, doublecomplex *q, int *ldq, double *rho, int *indxq, double *qstore, int *qptr, int *prmptr, int *perm, int *givptr, int *givcol, double *givnum, doublecomplex *work, double *rwork, int *iwork, int *info);

F77_RET_I zlaed8_f77(int *k, int *n, int *qsiz, doublecomplex *q, int *ldq, double *d, double *rho, int *cutpnt, double *z, double *dlamda, doublecomplex *q2, int *ldq2, double *w, int *indxp, int *indx, int *indxq, int *perm, int *givptr, int *givcol, double *givnum, int *info);

F77_RET_I zlaein_f77(logical *rightv, logical *noinit, int *n, doublecomplex *h, int *ldh, doublecomplex *w, doublecomplex *v, doublecomplex *b, int *ldb, double *rwork, double *eps3, double *smlnum, int *info);

F77_RET_I zlaesy_f77(doublecomplex *a, doublecomplex *b, doublecomplex *c, doublecomplex *rt1, doublecomplex *rt2, doublecomplex *evscal, doublecomplex *cs1, doublecomplex *sn1);

F77_RET_I zlaev2_f77(doublecomplex *a, doublecomplex *b, doublecomplex *c, double *rt1, double *rt2, double *cs1, doublecomplex *sn1);

F77_RET_I zlag2c_f77(int *m, int *n, doublecomplex *a, int *lda, singlecomplex *sa, int *ldsa, int *info);

F77_RET_I zlags2_f77(logical *upper, double *a1, doublecomplex *a2, double *a3, double *b1, doublecomplex *b2, double *b3, double *csu, doublecomplex *snu, double *csv, doublecomplex *snv, double *csq, doublecomplex *snq);

F77_RET_I zlagtm_f77(const char *trans, int *n, int *nrhs, double *alpha, doublecomplex *dl, doublecomplex *d, doublecomplex *du, doublecomplex *x, int *ldx, double *beta, doublecomplex *b, int *ldb);

F77_RET_I zlahef_f77(const char *uplo, int *n, int *nb, int *kb, doublecomplex *a, int *lda, int *ipiv, doublecomplex *w, int *ldw, int *info);

F77_RET_I zlahqr_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, doublecomplex *h, int *ldh, doublecomplex *w, int *iloz, int *ihiz, doublecomplex *z, int *ldz, int *info);

F77_RET_I zlahr2_f77(int *n, int *k, int *nb, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *t, int *ldt, doublecomplex *y, int *ldy);

F77_RET_I zlahrd_f77(int *n, int *k, int *nb, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *t, int *ldt, doublecomplex *y, int *ldy);

F77_RET_I zlaic1_f77(int *job, int *j, doublecomplex *x, double *sest, doublecomplex *w, doublecomplex *gamma, double *sestpr, doublecomplex *s, doublecomplex *c);

F77_RET_I zlals0_f77(int *icompq, int *nl, int *nr, int *sqre, int *nrhs, doublecomplex *b, int *ldb, doublecomplex *bx, int *ldbx, int *perm, int *givptr, int *givcol, int *ldgcol, double *givnum, int *ldgnum, double *poles, double *difl, double *difr, double *z, int *k, double *c, double *s, double *rwork, int *info);

F77_RET_I zlalsa_f77(int *icompq, int *smlsiz, int *n, int *nrhs, doublecomplex *b, int *ldb, doublecomplex *bx, int *ldbx, double *u, int *ldu, double *vt, int *k, double *difl, double *difr, double *z, double *poles, int *givptr, int *givcol, int *ldgcol, int *perm, double *givnum, double *c, double *s, double *rwork, int *iwork, int *info);

F77_RET_I zlalsd_f77(const char *uplo, int *smlsiz, int *n, int *nrhs, double *d, double *e, doublecomplex *b, int *ldb, double *rcond, int *rank, doublecomplex *work, double *rwork, int *iwork, int *info);

F77_RET_D zlangb_f77(const char *norm, int *n, int *kl, int *ku, doublecomplex *ab, int *ldab, double *work);

F77_RET_D zlange_f77(const char *norm, int *m, int *n, doublecomplex *a, int *lda, double *work);

F77_RET_D zlangt_f77(const char *norm, int *n, doublecomplex *dl, doublecomplex *d, doublecomplex *du);

F77_RET_D zlanhb_f77(const char *norm, const char *uplo, int *n, int *k, doublecomplex *ab, int *ldab, double *work);

F77_RET_D zlanhe_f77(const char *norm, const char *uplo, int *n, doublecomplex *a, int *lda, double *work);

F77_RET_D zlanhp_f77(const char *norm, const char *uplo, int *n, doublecomplex *ap, double *work);

F77_RET_D zlanhs_f77(const char *norm, int *n, doublecomplex *a, int *lda, double *work);

F77_RET_D zlanht_f77(const char *norm, int *n, double *d, doublecomplex *e);

F77_RET_D zlansb_f77(const char *norm, const char *uplo, int *n, int *k, doublecomplex *ab, int *ldab, double *work);

F77_RET_D zlansp_f77(const char *norm, const char *uplo, int *n, doublecomplex *ap, double *work);

F77_RET_D zlansy_f77(const char *norm, const char *uplo, int *n, doublecomplex *a, int *lda, double *work);

F77_RET_D zlantb_f77(const char *norm, const char *uplo, const char *diag, int *n, int *k, doublecomplex *ab, int *ldab, double *work);

F77_RET_D zlantp_f77(const char *norm, const char *uplo, const char *diag, int *n, doublecomplex *ap, double *work);

F77_RET_D zlantr_f77(const char *norm, const char *uplo, const char *diag, int *m, int *n, doublecomplex *a, int *lda, double *work);

F77_RET_I zlapll_f77(int *n, doublecomplex *x, int *incx, doublecomplex *y, int *incy, double *ssmin);

F77_RET_I zlapmt_f77(logical *forwrd, int *m, int *n, doublecomplex *x, int *ldx, int *k);

F77_RET_I zlaqgb_f77(int *m, int *n, int *kl, int *ku, doublecomplex *ab, int *ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, const char *equed);

F77_RET_I zlaqge_f77(int *m, int *n, doublecomplex *a, int *lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, const char *equed);

F77_RET_I zlaqhb_f77(const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, double *s, double *scond, double *amax, const char *equed);

F77_RET_I zlaqhe_f77(const char *uplo, int *n, doublecomplex *a, int *lda, double *s, double *scond, double *amax, const char *equed);

F77_RET_I zlaqhp_f77(const char *uplo, int *n, doublecomplex *ap, double *s, double *scond, double *amax, const char *equed);

F77_RET_I zlaqp2_f77(int *m, int *n, int *offset, doublecomplex *a, int *lda, int *jpvt, doublecomplex *tau, double *vn1, double *vn2, doublecomplex *work);

F77_RET_I zlaqps_f77(int *m, int *n, int *offset, int *nb, int *kb, doublecomplex *a, int *lda, int *jpvt, doublecomplex *tau, double *vn1, double *vn2, doublecomplex *auxv, doublecomplex *f, int *ldf);

F77_RET_I zlaqr0_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, doublecomplex *h, int *ldh, doublecomplex *w, int *iloz, int *ihiz, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, int *info);

F77_RET_I zlaqr1_f77(int *n, doublecomplex *h, int *ldh, doublecomplex *s1, doublecomplex *s2, doublecomplex *v);

F77_RET_I zlaqr2_f77(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, doublecomplex *h, int *ldh, int *iloz, int *ihiz, doublecomplex *z, int *ldz, int *ns, int *nd, doublecomplex *sh, doublecomplex *v, int *ldv, int *nh, doublecomplex *t, int *ldt, int *nv, doublecomplex *wv, int *ldwv, doublecomplex *work, int *lwork);

F77_RET_I zlaqr3_f77(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, doublecomplex *h, int *ldh, int *iloz, int *ihiz, doublecomplex *z, int *ldz, int *ns, int *nd, doublecomplex *sh, doublecomplex *v, int *ldv, int *nh, doublecomplex *t, int *ldt, int *nv, doublecomplex *wv, int *ldwv, doublecomplex *work, int *lwork);

F77_RET_I zlaqr4_f77(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, doublecomplex *h, int *ldh, doublecomplex *w, int *iloz, int *ihiz, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, int *info);

F77_RET_I zlaqr5_f77(logical *wantt, logical *wantz, int *kacc22, int *n, int *ktop, int *kbot, int *nshfts, doublecomplex *s, doublecomplex *h, int *ldh, int *iloz, int *ihiz, doublecomplex *z, int *ldz, doublecomplex *v, int *ldv, doublecomplex *u, int *ldu, int *nv, doublecomplex *wv, int *ldwv, int *nh, doublecomplex *wh, int *ldwh);

F77_RET_I zlaqsb_f77(const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, double *s, double *scond, double *amax, const char *equed);

F77_RET_I zlaqsp_f77(const char *uplo, int *n, doublecomplex *ap, double *s, double *scond, double *amax, const char *equed);

F77_RET_I zlaqsy_f77(const char *uplo, int *n, doublecomplex *a, int *lda, double *s, double *scond, double *amax, const char *equed);

F77_RET_I zlar1v_f77(int *n, int *b1, int *bn, double *lambda, double *d, double *l, double *ld, double *lld, double *pivmin, double *gaptol, doublecomplex *z, logical *wantnc, int *negcnt, double *ztz, double *mingma, int *r, int *isuppz, double *nrminv, double *resid, double *rqcorr, double *work);

F77_RET_I zlar2v_f77(int *n, doublecomplex *x, doublecomplex *y, doublecomplex *z, int *incx, double *c, doublecomplex *s, int *incc);

F77_RET_I zlarcm_f77(int *m, int *n, double *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *c, int *ldc, double *rwork);

F77_RET_I zlarf_f77(const char *side, int *m, int *n, doublecomplex *v, int *incv, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work);

F77_RET_I zlarfb_f77(const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n, int *k, doublecomplex *v, int *ldv, doublecomplex *t, int *ldt, doublecomplex *c, int *ldc, doublecomplex *work, int *ldwork);

F77_RET_I zlarfg_f77(int *n, doublecomplex *alpha, doublecomplex *x, int *incx, doublecomplex *tau);

F77_RET_I zlarft_f77(const char *direct, const char *storev, int *n, int *k, doublecomplex *v, int *ldv, doublecomplex *tau, doublecomplex *t, int *ldt);

F77_RET_I zlarfx_f77(const char *side, int *m, int *n, doublecomplex *v, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work);

F77_RET_I zlargv_f77(int *n, doublecomplex *x, int *incx, doublecomplex *y, int *incy, double *c, int *incc);

F77_RET_I zlarnv_f77(int *idist, int *iseed, int *n, doublecomplex *x);

F77_RET_I zlarrv_f77(int *n, double *vl, double *vu, double *d, double *l, double *pivmin, int *isplit, int *m, int *dol, int *dou, double *minrgp, double *rtol1, double *rtol2, double *w, double *werr, double *wgap, int *iblock, int *indexw, double *gers, doublecomplex *z, int *ldz, int *isuppz, double *work, int *iwork, int *info);

F77_RET_I zlartg_f77(doublecomplex *f, doublecomplex *g, double *cs, doublecomplex *sn, doublecomplex *r);

F77_RET_I zlartv_f77(int *n, doublecomplex *x, int *incx, doublecomplex *y, int *incy, double *c, doublecomplex *s, int *incc);

F77_RET_I zlarz_f77(const char *side, int *m, int *n, int *l, doublecomplex *v, int *incv, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work);

F77_RET_I zlarzb_f77(const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n, int *k, int *l, doublecomplex *v, int *ldv, doublecomplex *t, int *ldt, doublecomplex *c, int *ldc, doublecomplex *work, int *ldwork);

F77_RET_I zlarzt_f77(const char *direct, const char *storev, int *n, int *k, doublecomplex *v, int *ldv, doublecomplex *tau, doublecomplex *t, int *ldt);

F77_RET_I zlascl_f77(const char *type, int *kl, int *ku, double *cfrom, double *cto, int *m, int *n, doublecomplex *a, int *lda, int *info);

F77_RET_I zlaset_f77(const char *uplo, int *m, int *n, doublecomplex *alpha, doublecomplex *beta, doublecomplex *a, int *lda);

F77_RET_I zlasr_f77(const char *side, const char *pivot, const char *direct, int *m, int *n, double *c, double *s, doublecomplex *a, int *lda);

F77_RET_I zlassq_f77(int *n, doublecomplex *x, int *incx, double *scale, double *sumsq);

F77_RET_I zlaswp_f77(int *n, doublecomplex *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);

F77_RET_I zlasyf_f77(const char *uplo, int *n, int *nb, int *kb, doublecomplex *a, int *lda, int *ipiv, doublecomplex *w, int *ldw, int *info);

F77_RET_I zlatbs_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, int *kd, doublecomplex *ab, int *ldab, doublecomplex *x, double *scale, double *cnorm, int *info);

F77_RET_I zlatdf_f77(int *ijob, int *n, doublecomplex *z, int *ldz, doublecomplex *rhs, double *rdsum, double *rdscal, int *ipiv, int *jpiv);

F77_RET_I zlatps_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, doublecomplex *ap, doublecomplex *x, double *scale, double *cnorm, int *info);

F77_RET_I zlatrd_f77(const char *uplo, int *n, int *nb, doublecomplex *a, int *lda, double *e, doublecomplex *tau, doublecomplex *w, int *ldw);

F77_RET_I zlatrs_f77(const char *uplo, const char *trans, const char *diag, const char *normin, int *n, doublecomplex *a, int *lda, doublecomplex *x, double *scale, double *cnorm, int *info);

F77_RET_I zlatrz_f77(int *m, int *n, int *l, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work);

F77_RET_I zlatzm_f77(const char *side, int *m, int *n, doublecomplex *v, int *incv, doublecomplex *tau, doublecomplex *c1, doublecomplex *c2, int *ldc, doublecomplex *work);

F77_RET_I zlauu2_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *info);

F77_RET_I zlauum_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *info);

F77_RET_I zpbcon_f77(const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, double *anorm, double *rcond, doublecomplex *work, double *rwork, int *info);

F77_RET_I zpbequ_f77(const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, double *s, double *scond, double *amax, int *info);

F77_RET_I zpbrfs_f77(const char *uplo, int *n, int *kd, int *nrhs, doublecomplex *ab, int *ldab, doublecomplex *afb, int *ldafb, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zpbstf_f77(const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, int *info);

F77_RET_I zpbsv_f77(const char *uplo, int *n, int *kd, int *nrhs, doublecomplex *ab, int *ldab, doublecomplex *b, int *ldb, int *info);

F77_RET_I zpbsvx_f77(const char *fact, const char *uplo, int *n, int *kd, int *nrhs, doublecomplex *ab, int *ldab, doublecomplex *afb, int *ldafb, const char *equed, double *s, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zpbtf2_f77(const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, int *info);

F77_RET_I zpbtrf_f77(const char *uplo, int *n, int *kd, doublecomplex *ab, int *ldab, int *info);

F77_RET_I zpbtrs_f77(const char *uplo, int *n, int *kd, int *nrhs, doublecomplex *ab, int *ldab, doublecomplex *b, int *ldb, int *info);

F77_RET_I zpocon_f77(const char *uplo, int *n, doublecomplex *a, int *lda, double *anorm, double *rcond, doublecomplex *work, double *rwork, int *info);

F77_RET_I zpoequ_f77(int *n, doublecomplex *a, int *lda, double *s, double *scond, double *amax, int *info);

F77_RET_I zporfs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *af, int *ldaf, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zposv_f77(const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, int *info);

F77_RET_I zposvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *af, int *ldaf, const char *equed, double *s, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zpotf2_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *info);

F77_RET_I zpotrf_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *info);

F77_RET_I zpotri_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *info);

F77_RET_I zpotrs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, int *info);

F77_RET_I zppcon_f77(const char *uplo, int *n, doublecomplex *ap, double *anorm, double *rcond, doublecomplex *work, double *rwork, int *info);

F77_RET_I zppequ_f77(const char *uplo, int *n, doublecomplex *ap, double *s, double *scond, double *amax, int *info);

F77_RET_I zpprfs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *ap, doublecomplex *afp, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zppsv_f77(const char *uplo, int *n, int *nrhs, doublecomplex *ap, doublecomplex *b, int *ldb, int *info);

F77_RET_I zppsvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, doublecomplex *ap, doublecomplex *afp, const char *equed, double *s, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zpptrf_f77(const char *uplo, int *n, doublecomplex *ap, int *info);

F77_RET_I zpptri_f77(const char *uplo, int *n, doublecomplex *ap, int *info);

F77_RET_I zpptrs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *ap, doublecomplex *b, int *ldb, int *info);

F77_RET_I zptcon_f77(int *n, double *d, doublecomplex *e, double *anorm, double *rcond, double *rwork, int *info);

F77_RET_I zpteqr_f77(const char *compz, int *n, double *d, double *e, doublecomplex *z, int *ldz, double *work, int *info);

F77_RET_I zptrfs_f77(const char *uplo, int *n, int *nrhs, double *d, doublecomplex *e, double *df, doublecomplex *ef, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zptsv_f77(int *n, int *nrhs, double *d, doublecomplex *e, doublecomplex *b, int *ldb, int *info);

F77_RET_I zptsvx_f77(const char *fact, int *n, int *nrhs, double *d, doublecomplex *e, double *df, doublecomplex *ef, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zpttrf_f77(int *n, double *d, doublecomplex *e, int *info);

F77_RET_I zpttrs_f77(const char *uplo, int *n, int *nrhs, double *d, doublecomplex *e, doublecomplex *b, int *ldb, int *info);

F77_RET_I zptts2_f77(int *iuplo, int *n, int *nrhs, double *d, doublecomplex *e, doublecomplex *b, int *ldb);

F77_RET_I zrot_f77(int *n, doublecomplex *cx, int *incx, doublecomplex *cy, int *incy, double *c, doublecomplex *s);

F77_RET_I zspcon_f77(const char *uplo, int *n, doublecomplex *ap, int *ipiv, double *anorm, double *rcond, doublecomplex *work, int *info);

F77_RET_I zspmv_f77(const char *uplo, int *n, doublecomplex *alpha, doublecomplex *ap, doublecomplex *x, int *incx, doublecomplex *beta, doublecomplex *y, int *incy);

F77_RET_I zspr_f77(const char *uplo, int *n, doublecomplex *alpha, doublecomplex *x, int *incx, doublecomplex *ap);

F77_RET_I zsprfs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *ap, doublecomplex *afp, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zspsv_f77(const char *uplo, int *n, int *nrhs, doublecomplex *ap, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I zspsvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, doublecomplex *ap, doublecomplex *afp, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zsptrf_f77(const char *uplo, int *n, doublecomplex *ap, int *ipiv, int *info);

F77_RET_I zsptri_f77(const char *uplo, int *n, doublecomplex *ap, int *ipiv, doublecomplex *work, int *info);

F77_RET_I zsptrs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *ap, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I zstedc_f77(const char *compz, int *n, double *d, double *e, doublecomplex *z, int *ldz, doublecomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

F77_RET_I zstegr_f77(const char *jobz, const char *range, int *n, double *d, double *e, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, doublecomplex *z, int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I zstein_f77(int *n, double *d, double *e, int *m, double *w, int *iblock, int *isplit, doublecomplex *z, int *ldz, double *work, int *iwork, int *ifail, int *info);

F77_RET_I zstemr_f77(const char *jobz, const char *range, int *n, double *d, double *e, double *vl, double *vu, int *il, int *iu, int *m, double *w, doublecomplex *z, int *ldz, int *nzc, int *isuppz, logical *tryrac, double *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I zsteqr_f77(const char *compz, int *n, double *d, double *e, doublecomplex *z, int *ldz, double *work, int *info);

F77_RET_I zsycon_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *ipiv, double *anorm, double *rcond, doublecomplex *work, int *info);

F77_RET_I zsymv_f77(const char *uplo, int *n, doublecomplex *alpha, doublecomplex *a, int *lda, doublecomplex *x, int *incx, doublecomplex *beta, doublecomplex *y, int *incy);

F77_RET_I zsyr_f77(const char *uplo, int *n, doublecomplex *alpha, doublecomplex *x, int *incx, doublecomplex *a, int *lda);

F77_RET_I zsyrfs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *af, int *ldaf, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I zsysv_f77(const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *work, int *lwork, int *info);

F77_RET_I zsysvx_f77(const char *fact, const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *af, int *ldaf, int *ipiv, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *rcond, double *ferr, double *berr, doublecomplex *work, int *lwork, double *rwork, int *info);

F77_RET_I zsytf2_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *ipiv, int *info);

F77_RET_I zsytrf_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *ipiv, doublecomplex *work, int *lwork, int *info);

F77_RET_I zsytri_f77(const char *uplo, int *n, doublecomplex *a, int *lda, int *ipiv, doublecomplex *work, int *info);

F77_RET_I zsytrs_f77(const char *uplo, int *n, int *nrhs, doublecomplex *a, int *lda, int *ipiv, doublecomplex *b, int *ldb, int *info);

F77_RET_I ztbcon_f77(const char *norm, const char *uplo, const char *diag, int *n, int *kd, doublecomplex *ab, int *ldab, double *rcond, doublecomplex *work, double *rwork, int *info);

F77_RET_I ztbrfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *kd, int *nrhs, doublecomplex *ab, int *ldab, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I ztbtrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *kd, int *nrhs, doublecomplex *ab, int *ldab, doublecomplex *b, int *ldb, int *info);

F77_RET_I ztgevc_f77(const char *side, const char *howmny, logical *select, int *n, doublecomplex *s, int *lds, doublecomplex *p, int *ldp, doublecomplex *vl, int *ldvl, doublecomplex *vr, int *ldvr, int *mm, int *m, doublecomplex *work, double *rwork, int *info);

F77_RET_I ztgex2_f77(logical *wantq, logical *wantz, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *q, int *ldq, doublecomplex *z, int *ldz, int *j1, int *info);

F77_RET_I ztgexc_f77(logical *wantq, logical *wantz, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *q, int *ldq, doublecomplex *z, int *ldz, int *ifst, int *ilst, int *info);

F77_RET_I ztgsen_f77(int *ijob, logical *wantq, logical *wantz, logical *select, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *alpha, doublecomplex *beta, doublecomplex *q, int *ldq, doublecomplex *z, int *ldz, int *m, double *pl, double *pr, double *dif, doublecomplex *work, int *lwork, int *iwork, int *liwork, int *info);

F77_RET_I ztgsja_f77(const char *jobu, const char *jobv, const char *jobq, int *m, int *p, int *n, int *k, int *l, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, double *tola, double *tolb, double *alpha, double *beta, doublecomplex *u, int *ldu, doublecomplex *v, int *ldv, doublecomplex *q, int *ldq, doublecomplex *work, int *ncycle, int *info);

F77_RET_I ztgsna_f77(const char *job, const char *howmny, logical *select, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *vl, int *ldvl, doublecomplex *vr, int *ldvr, double *s, double *dif, int *mm, int *m, doublecomplex *work, int *lwork, int *iwork, int *info);

F77_RET_I ztgsy2_f77(const char *trans, int *ijob, int *m, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *c, int *ldc, doublecomplex *d, int *ldd, doublecomplex *e, int *lde, doublecomplex *f, int *ldf, double *scale, double *rdsum, double *rdscal, int *info);

F77_RET_I ztgsyl_f77(const char *trans, int *ijob, int *m, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *c, int *ldc, doublecomplex *d, int *ldd, doublecomplex *e, int *lde, doublecomplex *f, int *ldf, double *scale, double *dif, doublecomplex *work, int *lwork, int *iwork, int *info);

F77_RET_I ztpcon_f77(const char *norm, const char *uplo, const char *diag, int *n, doublecomplex *ap, double *rcond, doublecomplex *work, double *rwork, int *info);

F77_RET_I ztprfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, doublecomplex *ap, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I ztptri_f77(const char *uplo, const char *diag, int *n, doublecomplex *ap, int *info);

F77_RET_I ztptrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, doublecomplex *ap, doublecomplex *b, int *ldb, int *info);

F77_RET_I ztrcon_f77(const char *norm, const char *uplo, const char *diag, int *n, doublecomplex *a, int *lda, double *rcond, doublecomplex *work, double *rwork, int *info);

F77_RET_I ztrevc_f77(const char *side, const char *howmny, logical *select, int *n, doublecomplex *t, int *ldt, doublecomplex *vl, int *ldvl, doublecomplex *vr, int *ldvr, int *mm, int *m, doublecomplex *work, double *rwork, int *info);

F77_RET_I ztrexc_f77(const char *compq, int *n, doublecomplex *t, int *ldt, doublecomplex *q, int *ldq, int *ifst, int *ilst, int *info);

F77_RET_I ztrrfs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *x, int *ldx, double *ferr, double *berr, doublecomplex *work, double *rwork, int *info);

F77_RET_I ztrsen_f77(const char *job, const char *compq, logical *select, int *n, doublecomplex *t, int *ldt, doublecomplex *q, int *ldq, doublecomplex *w, int *m, double *s, double *sep, doublecomplex *work, int *lwork, int *info);

F77_RET_I ztrsna_f77(const char *job, const char *howmny, logical *select, int *n, doublecomplex *t, int *ldt, doublecomplex *vl, int *ldvl, doublecomplex *vr, int *ldvr, double *s, double *sep, int *mm, int *m, doublecomplex *work, int *ldwork, double *rwork, int *info);

F77_RET_I ztrsyl_f77(const char *trana, const char *tranb, int *isgn, int *m, int *n, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, doublecomplex *c, int *ldc, double *scale, int *info);

F77_RET_I ztrti2_f77(const char *uplo, const char *diag, int *n, doublecomplex *a, int *lda, int *info);

F77_RET_I ztrtri_f77(const char *uplo, const char *diag, int *n, doublecomplex *a, int *lda, int *info);

F77_RET_I ztrtrs_f77(const char *uplo, const char *trans, const char *diag, int *n, int *nrhs, doublecomplex *a, int *lda, doublecomplex *b, int *ldb, int *info);

F77_RET_I ztzrqf_f77(int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, int *info);

F77_RET_I ztzrzf_f77(int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zung2l_f77(int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *info);

F77_RET_I zung2r_f77(int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *info);

F77_RET_I zungbr_f77(const char *vect, int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zunghr_f77(int *n, int *ilo, int *ihi, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zungl2_f77(int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *info);

F77_RET_I zunglq_f77(int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zungql_f77(int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zungqr_f77(int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zungr2_f77(int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *info);

F77_RET_I zungrq_f77(int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zungtr_f77(const char *uplo, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *work, int *lwork, int *info);

F77_RET_I zunm2l_f77(const char *side, const char *trans, int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *info);

F77_RET_I zunm2r_f77(const char *side, const char *trans, int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *info);

F77_RET_I zunmbr_f77(const char *vect, const char *side, const char *trans, int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *lwork, int *info);

F77_RET_I zunmhr_f77(const char *side, const char *trans, int *m, int *n, int *ilo, int *ihi, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *lwork, int *info);

F77_RET_I zunml2_f77(const char *side, const char *trans, int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *info);

F77_RET_I zunmlq_f77(const char *side, const char *trans, int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *lwork, int *info);

F77_RET_I zunmql_f77(const char *side, const char *trans, int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *lwork, int *info);

F77_RET_I zunmqr_f77(const char *side, const char *trans, int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *lwork, int *info);

F77_RET_I zunmr2_f77(const char *side, const char *trans, int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *info);

F77_RET_I zunmr3_f77(const char *side, const char *trans, int *m, int *n, int *k, int *l, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *info);

F77_RET_I zunmrq_f77(const char *side, const char *trans, int *m, int *n, int *k, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *lwork, int *info);

F77_RET_I zunmrz_f77(const char *side, const char *trans, int *m, int *n, int *k, int *l, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *lwork, int *info);

F77_RET_I zunmtr_f77(const char *side, const char *uplo, const char *trans, int *m, int *n, doublecomplex *a, int *lda, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *lwork, int *info);

F77_RET_I zupgtr_f77(const char *uplo, int *n, doublecomplex *ap, doublecomplex *tau, doublecomplex *q, int *ldq, doublecomplex *work, int *info);

F77_RET_I zupmtr_f77(const char *side, const char *uplo, const char *trans, int *m, int *n, doublecomplex *ap, doublecomplex *tau, doublecomplex *c, int *ldc, doublecomplex *work, int *info);

#ifdef __cplusplus
}
#endif

#endif
