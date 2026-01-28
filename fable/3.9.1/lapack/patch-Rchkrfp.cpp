--- Rchkrfp.cpp_	2026-01-28 09:31:30.439494363 +0900
+++ Rchkrfp.cpp	2026-01-28 09:31:30.445494494 +0900
@@ -43,8 +43,10 @@
 #include <mplapack_matgen.h>
 #include <mplapack_lin.h>
 
-void program_dchkrfp(int argc, char const *argv[]) {
-    common cmn(argc, argv);
+#include <memory>
+
+void Rchkrfp(void) {
+    common cmn;
     common_read read(cmn);
     common_write write(cmn);
     //
@@ -53,9 +55,9 @@
     static const char *format_9997 = "(' Total time used = ',f12.2,' seconds',/)";
     static const char *format_9996 = "(' !! Invalid input value: ',a4,'=',i6,'; must be >=',i6)";
     static const char *format_9995 = "(' !! Invalid input value: ',a4,'=',i6,'; must be <=',i6)";
-    static const char *format_9994 = "(/,' Tests of the DOUBLE PRECISION LAPACK RFP routines ',/,"
-                                     "' LAPACK VERSION ',i1,'.',i1,'.',i1,/,/,"
-                                     "' The following parameter values will be used:')";
+    static const char *format_9994 = "(' Tests of the Multiple precision version of LAPACK RFP routines',i1,'.',i1,'.',i1,/, "
+                                     "' Based on original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, "
+                                     "'The following parameter values will be used:')";
     static const char *format_9993 = "(4x,a4,':  ',10i6,/,11x,10i6)";
     static const char *format_9992 = "(/,' Routines pass computational tests if test ratio is ','less than',"
                                      "f8.2,/)";
@@ -68,16 +70,19 @@
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
@@ -236,23 +241,40 @@
     // Test the routines: dpftrf, dpftri, dpftrs (as in Rdrvpo).
     // This also tests the routines: dtfsm, dtftri, dtfttr, dtrttf.
     //
-    REAL worka[nmax * nmax];
-    REAL workasav[nmax * nmax];
-    REAL workafac[nmax * nmax];
-    REAL workainv[nmax * nmax];
-    REAL workb[nmax * maxrhs];
-    REAL workbsav[nmax * maxrhs];
-    REAL workxact[nmax * maxrhs];
-    REAL workx[nmax * maxrhs];
-    REAL workarf[(nmax * (nmax + 1)) / 2];
-    REAL workarfinv[(nmax * (nmax + 1)) / 2];
-    REAL d_work_dlatms[3 * nmax];
-    REAL d_work_dpot01[nmax];
-    REAL d_temp_dpot02[nmax * maxrhs];
-    REAL d_temp_dpot03[nmax * nmax];
-    REAL d_work_dlansy[nmax];
-    REAL d_work_dpot02[nmax];
-    REAL d_work_dpot03[nmax];
+    auto worka_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * nmax));
+    REAL *worka = worka_storage.get();
+    auto workasav_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * nmax));
+    REAL *workasav = workasav_storage.get();
+    auto workafac_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * nmax));
+    REAL *workafac = workafac_storage.get();
+    auto workainv_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * nmax));
+    REAL *workainv = workainv_storage.get();
+    auto workb_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs));
+    REAL *workb = workb_storage.get();
+    auto workbsav_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs));
+    REAL *workbsav = workbsav_storage.get();
+    auto workxact_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs));
+    REAL *workxact = workxact_storage.get();
+    auto workx_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs));
+    REAL *workx = workx_storage.get();
+    auto workarf_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, (nmax * (nmax + 1)) / 2));
+    REAL *workarf = workarf_storage.get();
+    auto workarfinv_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, (nmax * (nmax + 1)) / 2));
+    REAL *workarfinv = workarfinv_storage.get();
+    auto d_work_dlatms_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, 3 * nmax));
+    REAL *d_work_dlatms = d_work_dlatms_storage.get();
+    auto d_work_dpot01_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax));
+    REAL *d_work_dpot01 = d_work_dpot01_storage.get();
+    auto d_temp_dpot02_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs));
+    REAL *d_temp_dpot02 = d_temp_dpot02_storage.get();
+    auto d_temp_dpot03_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * nmax));
+    REAL *d_temp_dpot03 = d_temp_dpot03_storage.get();
+    auto d_work_dlansy_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax));
+    REAL *d_work_dlansy = d_work_dlansy_storage.get();
+    auto d_work_dpot02_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax));
+    REAL *d_work_dpot02 = d_work_dpot02_storage.get();
+    auto d_work_dpot03_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax));
+    REAL *d_work_dpot03 = d_work_dpot03_storage.get();
     Rdrvrfp(nout, nn, nval, nns, nsval, nnt, ntval, thresh, worka, workasav, workafac, workainv, workb, workbsav, workxact, workx, workarf, workarfinv, d_work_dlatms, d_work_dpot01, d_temp_dpot02, d_temp_dpot03, d_work_dlansy, d_work_dpot02, d_work_dpot03);
     //
     // Test the routine: dlansf
@@ -262,7 +284,8 @@
     // Test the conversion routines:
     // dtfttp, dtpttf, dtfttr, dtrttf, dtrttp and dtpttr.
     //
-    REAL workap[(nmax * (nmax + 1)) / 2];
+    auto workap_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, (nmax * (nmax + 1)) / 2));
+    REAL* workap = workap_storage.get();
     Rdrvrf2(nout, nn, nval, worka, nmax, workarf, workap, workasav);
     //
     // Test the routine: dtfsm
@@ -282,4 +305,4 @@
     //
 }
 
-int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dchkrfp); }
+int main(int argc, char const *argv[]) { Rchkrfp(); }
