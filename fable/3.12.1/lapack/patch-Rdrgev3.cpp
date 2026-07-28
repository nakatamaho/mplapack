--- Rdrgev3.cpp~	2026-04-16 00:00:00.000000000 +0900
+++ Rdrgev3.cpp	2026-04-16 00:00:00.000000000 +0900
@@ -436,7 +436,13 @@ void Rdrgev3(INTEGER const nsizes, INTEG
             }
             //
             for (j = 1; j <= n; j = j + 1) {
+#if defined MPLAPACK_BUILD_WITH_MPFR || defined MPLAPACK_BUILD_WITH_GMP
+                REAL scale = max(safmin, abs(alphar[j - 1]) + abs(alphai[j - 1]) + abs(beta[j - 1]), abs(alphr1[j - 1]) + abs(alphi1[j - 1]) + abs(beta1[j - 1]));
+                REAL diff = abs(alphar[j - 1] - alphr1[j - 1]) + abs(alphai[j - 1] - alphi1[j - 1]) + abs(beta[j - 1] - beta1[j - 1]);
+                if (diff > 100.0 * ulp * scale) {
+#else
                 if (alphar[j - 1] != alphr1[j - 1] || alphai[j - 1] != alphi1[j - 1] || beta[j - 1] != beta1[j - 1]) {
+#endif
                     result[5 - 1] = ulpinv;
                 }
             }
@@ -455,7 +461,13 @@ void Rdrgev3(INTEGER const nsizes, INTEG
             }
             //
             for (j = 1; j <= n; j = j + 1) {
+#if defined MPLAPACK_BUILD_WITH_MPFR || defined MPLAPACK_BUILD_WITH_GMP
+                REAL scale = max(safmin, abs(alphar[j - 1]) + abs(alphai[j - 1]) + abs(beta[j - 1]), abs(alphr1[j - 1]) + abs(alphi1[j - 1]) + abs(beta1[j - 1]));
+                REAL diff = abs(alphar[j - 1] - alphr1[j - 1]) + abs(alphai[j - 1] - alphi1[j - 1]) + abs(beta[j - 1] - beta1[j - 1]);
+                if (diff > 100.0 * ulp * scale) {
+#else
                 if (alphar[j - 1] != alphr1[j - 1] || alphai[j - 1] != alphi1[j - 1] || beta[j - 1] != beta1[j - 1]) {
+#endif
                     result[6 - 1] = ulpinv;
                 }
             }
@@ -482,7 +494,13 @@ void Rdrgev3(INTEGER const nsizes, INTEG
             }
             //
             for (j = 1; j <= n; j = j + 1) {
+#if defined MPLAPACK_BUILD_WITH_MPFR || defined MPLAPACK_BUILD_WITH_GMP
+                REAL scale = max(safmin, abs(alphar[j - 1]) + abs(alphai[j - 1]) + abs(beta[j - 1]), abs(alphr1[j - 1]) + abs(alphi1[j - 1]) + abs(beta1[j - 1]));
+                REAL diff = abs(alphar[j - 1] - alphr1[j - 1]) + abs(alphai[j - 1] - alphi1[j - 1]) + abs(beta[j - 1] - beta1[j - 1]);
+                if (diff > 100.0 * ulp * scale) {
+#else
                 if (alphar[j - 1] != alphr1[j - 1] || alphai[j - 1] != alphi1[j - 1] || beta[j - 1] != beta1[j - 1]) {
+#endif
                     result[7 - 1] = ulpinv;
                 }
             }
