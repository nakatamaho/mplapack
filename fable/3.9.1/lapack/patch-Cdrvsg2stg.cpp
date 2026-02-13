--- Cdrvsg2stg.cpp
+++ Cdrvsg2stg.cpp
@@ -637,7 +637,7 @@ void Cdrvsg2stg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *d
                         }
                     }
                     //
-                    Chpgvx(ibtype, "V", "A", uplo.elems, n, ap, bp, vl, vu, il, iu, abstol, m, d, z, ldz, work, rwork, &iwork[(n + 1) - 1], iwork, info);
+                    Chpgvx(ibtype, "V", "A", uplo.elems, n, ap, bp, vl, vu, il, iu, abstol, m, d, z, ldz, work, rwork, &iwork[(n + 1) - 1], iwork, iinfo);
                     if (iinfo != 0) {
                         write(nounit, format_9999), "Chpgvx(V,A" + uplo + ")", iinfo, n, jtype, ioldsd;
                         info = abs(iinfo);
@@ -679,7 +679,7 @@ void Cdrvsg2stg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *d
                     //
                     vl = zero;
                     vu = anorm;
-                    Chpgvx(ibtype, "V", "V", uplo.elems, n, ap, bp, vl, vu, il, iu, abstol, m, d, z, ldz, work, rwork, &iwork[(n + 1) - 1], iwork, info);
+                    Chpgvx(ibtype, "V", "V", uplo.elems, n, ap, bp, vl, vu, il, iu, abstol, m, d, z, ldz, work, rwork, &iwork[(n + 1) - 1], iwork, iinfo);
                     if (iinfo != 0) {
                         write(nounit, format_9999), "Chpgvx(V,V" + uplo + ")", iinfo, n, jtype, ioldsd;
                         info = abs(iinfo);
@@ -719,7 +719,7 @@ void Cdrvsg2stg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *d
                         }
                     }
                     //
-                    Chpgvx(ibtype, "V", "I", uplo.elems, n, ap, bp, vl, vu, il, iu, abstol, m, d, z, ldz, work, rwork, &iwork[(n + 1) - 1], iwork, info);
+                    Chpgvx(ibtype, "V", "I", uplo.elems, n, ap, bp, vl, vu, il, iu, abstol, m, d, z, ldz, work, rwork, &iwork[(n + 1) - 1], iwork, iinfo);
                     if (iinfo != 0) {
                         write(nounit, format_9999), "Chpgvx(V,I" + uplo + ")", iinfo, n, jtype, ioldsd;
                         info = abs(iinfo);
