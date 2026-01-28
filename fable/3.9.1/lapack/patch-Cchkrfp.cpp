--- Cchkrfp.cpp_	2026-01-28 09:31:30.411493758 +0900
+++ Cchkrfp.cpp	2026-01-28 09:31:30.417493887 +0900
@@ -43,8 +43,10 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void program_zchkrfp(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+#include <memory>
+
+void Cchkrfp(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     //
@@ -53,8 +55,9 @@
     static const char *format_9997 = "(' Total time used = ',f12.2,' seconds',/)";
     static const char *format_9996 = "(' !! Invalid input value: ',a4,'=',i6,'; must be >=',i6)";
     static const char *format_9995 = "(' !! Invalid input value: ',a4,'=',i6,'; must be <=',i6)";
-    static const char *format_9994 = "(/,' Tests of the COMPLEX*16 LAPACK RFP routines ',/,' LAPACK VERSION ',"
-                                     "i1,'.',i1,'.',i1,/,/,' The following parameter values will be used:')";
+    static const char *format_9994 = "(' Tests of the Multiple precision version of LAPACK MPLAPACK VERSION ',i1,'.',i1,'.',i1,/, "
+                                     "' Based on original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, "
+                                     "'The following parameter values will be used:')";
     static const char *format_9993 = "(4x,a4,':  ',10i6,/,11x,10i6)";
     static const char *format_9992 = "(/,' Routines pass computational tests if test ratio is ','less than',"
                                      "f8.2,/)";
@@ -67,16 +70,19 @@
     // Read a dummy line.
     //
     const INTEGER nin = 5;
+    const INTEGER nout = 6;
     read(nin, star);
     //
     // Report LAPACK version tag (e.g. LAPACK-3.2.0)
     //
-    INTEGER vers_major = 0;
-    INTEGER vers_minor = 0;
-    INTEGER vers_patch = 0;
-    ilaver(vers_major, vers_minor, vers_patch);
-    const INTEGER nout = 6;
-    write(nout, format_9994), vers_major, vers_minor, vers_patch;
+    INTEGER mplapack_vers_major = 0;
+    INTEGER mplapack_vers_minor = 0;
+    INTEGER mplapack_vers_patch = 0;
+    INTEGER lapack_vers_major = 0;
+    INTEGER lapack_vers_minor = 0;
+    INTEGER lapack_vers_patch = 0;
+    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
+    write(nout, format_9994), mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
     //
     // Read the values of N
     //
@@ -235,24 +241,42 @@
     // Test the routines: zpftrf, zpftri, zpftrs (as in Cdrvpo).
     // This also tests the routines: ztfsm, ztftri, ztfttr, ztrttf.
     //
-    COMPLEX worka[nmax * nmax];
-    COMPLEX workasav[nmax * nmax];
-    COMPLEX workafac[nmax * nmax];
-    COMPLEX workainv[nmax * nmax];
-    COMPLEX workb[nmax * maxrhs];
-    COMPLEX workbsav[nmax * maxrhs];
-    COMPLEX workxact[nmax * maxrhs];
-    COMPLEX workx[nmax * maxrhs];
-    COMPLEX workarf[(nmax * (nmax + 1)) / 2];
-    COMPLEX workarfinv[(nmax * (nmax + 1)) / 2];
-    COMPLEX z_work_zlatms[3 * nmax];
-    COMPLEX z_work_zpot02[nmax * maxrhs];
-    COMPLEX z_work_zpot03[nmax * nmax];
-    REAL d_work_zlatms[nmax];
-    REAL d_work_zlanhe[nmax];
-    REAL d_work_zpot01[nmax];
-    REAL d_work_zpot02[nmax];
-    REAL d_work_zpot03[nmax];
+    auto worka_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * nmax));
+    COMPLEX *worka = worka_storage.get();
+    auto workasav_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * nmax));
+    COMPLEX *workasav = workasav_storage.get();
+    auto workafac_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * nmax));
+    COMPLEX *workafac = workafac_storage.get();
+    auto workainv_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * nmax));
+    COMPLEX *workainv = workainv_storage.get();
+    auto workb_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * maxrhs));
+    COMPLEX *workb = workb_storage.get();
+    auto workbsav_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * maxrhs));
+    COMPLEX *workbsav = workbsav_storage.get();
+    auto workxact_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * maxrhs));
+    COMPLEX *workxact = workxact_storage.get();
+    auto workx_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * maxrhs));
+    COMPLEX *workx = workx_storage.get();
+    auto workarf_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, (nmax * (nmax + 1)) / 2));
+    COMPLEX *workarf = workarf_storage.get();
+    auto workarfinv_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, (nmax * (nmax + 1)) / 2));
+    COMPLEX *workarfinv = workarfinv_storage.get();
+    auto z_work_zlatms_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, 3 * nmax));
+    COMPLEX *z_work_zlatms = z_work_zlatms_storage.get();
+    auto z_work_zpot02_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * maxrhs));
+    COMPLEX *z_work_zpot02 = z_work_zpot02_storage.get();
+    auto z_work_zpot03_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * nmax));
+    COMPLEX *z_work_zpot03 = z_work_zpot03_storage.get();
+    auto d_work_zlatms_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax));
+    REAL *d_work_zlatms = d_work_zlatms_storage.get();
+    auto d_work_zlanhe_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax));
+    REAL *d_work_zlanhe = d_work_zlanhe_storage.get();
+    auto d_work_zpot01_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax));
+    REAL *d_work_zpot01 = d_work_zpot01_storage.get();
+    auto d_work_zpot02_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax));
+    REAL *d_work_zpot02 = d_work_zpot02_storage.get();
+    auto d_work_zpot03_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax));
+    REAL *d_work_zpot03 = d_work_zpot03_storage.get();
     Cdrvrfp(nout, nn, nval, nns, nsval, nnt, ntval, thresh, worka, workasav, workafac, workainv, workb, workbsav, workxact, workx, workarf, workarfinv, z_work_zlatms, z_work_zpot02, z_work_zpot03, d_work_zlatms, d_work_zlanhe, d_work_zpot01, d_work_zpot02, d_work_zpot03);
     //
     // Test the routine: zlanhf
@@ -262,7 +286,8 @@
     // Test the conversion routines:
     // zhfttp, ztpthf, ztfttr, ztrttf, ztrttp and ztpttr.
     //
-    COMPLEX workap[(nmax * (nmax + 1)) / 2];
+    auto workap_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, (nmax * (nmax + 1)) / 2));
+    COMPLEX *workap = workap_storage.get();
     Cdrvrf2(nout, nn, nval, worka, nmax, workarf, workap, workasav);
     //
     // Test the routine: ztfsm
@@ -282,4 +307,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_zchkrfp); }
+int main(int argc, char const *argv[]) { Cchkrfp(); }
