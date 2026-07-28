--- Rget32.cpp~	2026-04-16 00:00:00.000000000 +0900
+++ Rget32.cpp	2026-04-16 00:00:00.000000000 +0900
@@ -273,10 +273,28 @@ void Rget32(REAL &rmax, INTEGER &lmax, I
                                             }
                                             tnrm = abs(tr[0]) + abs(tr[(2 - 1)]) + abs(tr[(2 - 1) * ldtr]) + abs(tr[(2 - 1) + (2 - 1) * ldtr]) + abs(tl[0]) + abs(tl[(2 - 1)]) + abs(tl[(2 - 1) * ldtl]) + abs(tl[(2 - 1) + (2 - 1) * ldtl]);
                                             xnrm = max(abs(x[0]) + abs(x[(2 - 1) * ldx]), abs(x[(2 - 1)]) + abs(x[(2 - 1) + (2 - 1) * ldx]));
+#if defined MPLAPACK_BUILD_WITH_QD
+                                            REAL bsum = abs(b[0]) + abs(b[(2 - 1)]) + abs(b[(2 - 1) * ldb]) + abs(b[(2 - 1) + (2 - 1) * ldb]);
+                                            REAL cscale = max(tnrm, smlnum);
+                                            REAL xscale = max(xnrm, abs(scale) * bsum / cscale, smlnum);
+                                            REAL rscale = cscale * xscale;
+                                            REAL cinv = one / cscale;
+                                            REAL xinv = one / xscale;
+                                            REAL rinv = one / rscale;
+                                            REAL res1_qd = rscale * abs(((tl[0] * cinv + sgn * tr[0] * cinv)) * (x[0] * xinv) + (sgn * tr[(2 - 1)] * cinv) * (x[(2 - 1) * ldx] * xinv) + (tl[(2 - 1) * ldtl] * cinv) * (x[(2 - 1)] * xinv) - (scale * b[0] * rinv));
+                                            REAL res2_qd = rscale * abs((tl[0] * cinv) * (x[(2 - 1) * ldx] * xinv) + (sgn * tr[(2 - 1) * ldtr] * cinv) * (x[0] * xinv) + (sgn * tr[(2 - 1) + (2 - 1) * ldtr] * cinv) * (x[(2 - 1) * ldx] * xinv) + (tl[(2 - 1) * ldtl] * cinv) * (x[(2 - 1) + (2 - 1) * ldx] * xinv) - (scale * b[(2 - 1) * ldb] * rinv));
+                                            REAL res3_qd = rscale * abs((tl[(2 - 1)] * cinv) * (x[0] * xinv) + (sgn * tr[0] * cinv) * (x[(2 - 1)] * xinv) + (sgn * tr[(2 - 1)] * cinv) * (x[(2 - 1) + (2 - 1) * ldx] * xinv) + (tl[(2 - 1) + (2 - 1) * ldtl] * cinv) * (x[(2 - 1)] * xinv) - (scale * b[(2 - 1)] * rinv));
+                                            REAL res4_qd = rscale * abs(((tl[(2 - 1) + (2 - 1) * ldtl] * cinv + sgn * tr[(2 - 1) + (2 - 1) * ldtr] * cinv)) * (x[(2 - 1) + (2 - 1) * ldx] * xinv) + (sgn * tr[(2 - 1) * ldtr] * cinv) * (x[(2 - 1)] * xinv) + (tl[(2 - 1)] * cinv) * (x[(2 - 1) * ldx] * xinv) - (scale * b[(2 - 1) + (2 - 1) * ldb] * rinv));
+                                            res = res1_qd;
+                                            res += res2_qd;
+                                            res += res3_qd;
+                                            res += res4_qd;
+#else
                                             res = abs(((tl[0] + sgn * tr[0])) * (x[0]) + (sgn * tr[(2 - 1)]) * (x[(2 - 1) * ldx]) + (tl[(2 - 1) * ldtl]) * (x[(2 - 1)]) - (scale * b[0]));
                                             res += abs((tl[0]) * (x[(2 - 1) * ldx]) + (sgn * tr[(2 - 1) * ldtr]) * (x[0]) + (sgn * tr[(2 - 1) + (2 - 1) * ldtr]) * (x[(2 - 1) * ldx]) + (tl[(2 - 1) * ldtl]) * (x[(2 - 1) + (2 - 1) * ldx]) - (scale * b[(2 - 1) * ldb]));
                                             res += abs((tl[(2 - 1)]) * (x[0]) + (sgn * tr[0]) * (x[(2 - 1)]) + (sgn * tr[(2 - 1)]) * (x[(2 - 1) + (2 - 1) * ldx]) + (tl[(2 - 1) + (2 - 1) * ldtl]) * (x[(2 - 1)]) - (scale * b[(2 - 1)]));
                                             res += abs(((tl[(2 - 1) + (2 - 1) * ldtl] + sgn * tr[(2 - 1) + (2 - 1) * ldtr])) * (x[(2 - 1) + (2 - 1) * ldx]) + (sgn * tr[(2 - 1) * ldtr]) * (x[(2 - 1)]) + (tl[(2 - 1)]) * (x[(2 - 1) * ldx]) - (scale * b[(2 - 1) + (2 - 1) * ldb]));
+#endif
                                             den = max(smlnum, smlnum * xnrm, (tnrm * eps) * xnrm);
                                             res = res / den;
                                             if (scale > one) {
