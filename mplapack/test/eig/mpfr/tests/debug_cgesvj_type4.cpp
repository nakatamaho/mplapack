// debug_cgesvj_type4_v2.cpp
// Revised: uses Clatms (or manual unitary transform) to generate
// dense type 4 matrix matching what Cdrvbd actually tests.
//
// Build example:
//   g++ -O2 -I${MPLAPACK}/include \
//       -D___MPLAPACK_BUILD_WITH_MPFR___ \
//       debug_cgesvj_type4_v2.cpp \
//       -L${MPLAPACK}/lib -lmplapack_mpfr -lmpblas_mpfr \
//       -lmpc -lmpfr -lgmp -o debug_cgesvj_type4_v2
//
// If Clatms is not easily callable, this program uses an alternative:
// apply Householder reflections to a diagonal matrix to make it dense.

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>
#include <cstdio>
#include <cstdlib>

// ---------- helpers ----------
void _printnum(mpreal a) {
    mpfr_printf("%10.8Re", mpfr_ptr(a));
}

void _printnum_full(mpreal a) {
    mpfr_printf("%40.34Re", mpfr_ptr(a));
}

// Compute Frobenius norm using Rlamch-aware scaling to avoid underflow
// (mimics Rlange("F",...))
mpreal safe_frobenius_norm_complex(int M, int N, mpcomplex *A, int LDA) {
    mpreal scale = 0.0;
    mpreal ssq   = 1.0;
    mpreal zero  = 0.0;
    // Use scaled sum of squares: ||A||_F = scale * sqrt(ssq)
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            mpreal re = abs(A[i + j * LDA].real());
            mpreal im = abs(A[i + j * LDA].imag());
            // Process real part
            if (re != zero) {
                if (scale < re) {
                    ssq = 1.0 + ssq * (scale / re) * (scale / re);
                    scale = re;
                } else {
                    ssq += (re / scale) * (re / scale);
                }
            }
            // Process imaginary part
            if (im != zero) {
                if (scale < im) {
                    ssq = 1.0 + ssq * (scale / im) * (scale / im);
                    scale = im;
                } else {
                    ssq += (im / scale) * (im / scale);
                }
            }
        }
    }
    return scale * sqrt(ssq);
}

// Apply a random Householder reflection from left: A <- (I - 2 v v^H) * A
// v is a random unit vector of length M
void apply_random_left_householder(int M, int N, mpcomplex *A, int LDA,
                                    unsigned long seed) {
    mpreal one = 1.0, zero = 0.0, two = 2.0;
    mpcomplex *v = new mpcomplex[M];
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);

    // Generate random vector
    mpreal norm_sq = zero;
    for (int i = 0; i < M; i++) {
        mpfr_t re_t, im_t;
        mpfr_init2(re_t, mpfr_get_default_prec());
        mpfr_init2(im_t, mpfr_get_default_prec());
        mpfr_urandomb(re_t, state);
        mpfr_urandomb(im_t, state);
        mpreal re_val(re_t), im_val(im_t);
        re_val = re_val - 0.5;
        im_val = im_val - 0.5;
        v[i] = mpcomplex(re_val, im_val);
        norm_sq += re_val * re_val + im_val * im_val;
        mpfr_clear(re_t);
        mpfr_clear(im_t);
    }

    // Normalize v
    mpreal norm_v = sqrt(norm_sq);
    for (int i = 0; i < M; i++) {
        v[i] = v[i] / mpcomplex(norm_v, zero);
    }

    // Apply: A <- A - 2 * v * (v^H * A)
    // First compute w = v^H * A (1 x N)
    mpcomplex *w = new mpcomplex[N];
    for (int j = 0; j < N; j++) {
        w[j] = mpcomplex(zero, zero);
        for (int i = 0; i < M; i++) {
            w[j] += conj(v[i]) * A[i + j * LDA];
        }
    }
    // Then A <- A - 2 * v * w^T
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            A[i + j * LDA] -= two * v[i] * w[j];
        }
    }

    delete[] v;
    delete[] w;
    gmp_randclear(state);
}

// Apply a random Householder reflection from right: A <- A * (I - 2 v v^H)
void apply_random_right_householder(int M, int N, mpcomplex *A, int LDA,
                                     unsigned long seed) {
    mpreal one = 1.0, zero = 0.0, two = 2.0;
    mpcomplex *v = new mpcomplex[N];
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);

    // Generate random vector
    mpreal norm_sq = zero;
    for (int i = 0; i < N; i++) {
        mpfr_t re_t, im_t;
        mpfr_init2(re_t, mpfr_get_default_prec());
        mpfr_init2(im_t, mpfr_get_default_prec());
        mpfr_urandomb(re_t, state);
        mpfr_urandomb(im_t, state);
        mpreal re_val(re_t), im_val(im_t);
        re_val = re_val - 0.5;
        im_val = im_val - 0.5;
        v[i] = mpcomplex(re_val, im_val);
        norm_sq += re_val * re_val + im_val * im_val;
        mpfr_clear(re_t);
        mpfr_clear(im_t);
    }
    mpreal norm_v = sqrt(norm_sq);
    for (int i = 0; i < N; i++) {
        v[i] = v[i] / mpcomplex(norm_v, zero);
    }

    // Apply: A <- A - 2 * (A * v) * v^H
    // First compute w = A * v (M x 1)
    mpcomplex *w = new mpcomplex[M];
    for (int i = 0; i < M; i++) {
        w[i] = mpcomplex(zero, zero);
        for (int j = 0; j < N; j++) {
            w[i] += A[i + j * LDA] * v[j];
        }
    }
    // Then A <- A - 2 * w * v^H
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            A[i + j * LDA] -= two * w[i] * conj(v[j]);
        }
    }

    delete[] v;
    delete[] w;
    gmp_randclear(state);
}

// ---------- main ----------
int main() {
    printf("==================================================\n");
    printf("  Cgesvj debug test v2: type 4 (dense matrix)\n");
    printf("==================================================\n\n");

    // ---- 1. MPFR profile ----
    printf("--- MPFR profile ---\n");
    printf("  mpfr_get_default_prec()  = %ld\n", (long)mpfr_get_default_prec());
    printf("  mpfr_get_emin()          = %ld\n", (long)mpfr_get_emin());
    printf("  mpfr_get_emax()          = %ld\n", (long)mpfr_get_emax());

    mpreal eps    = Rlamch_mpfr("E");
    mpreal sfmin  = Rlamch_mpfr("S");
    mpreal ovfl   = Rlamch_mpfr("O");
    mpreal unfl   = Rlamch_mpfr("U");
    mpreal one    = 1.0;
    mpreal zero   = 0.0;

    printf("  Rlamch('E') = "); _printnum_full(eps);   printf("\n");
    printf("  Rlamch('S') = "); _printnum_full(sfmin); printf("\n");
    printf("  Rlamch('O') = "); _printnum_full(ovfl);  printf("\n");

    mpreal rootsfmin = sqrt(sfmin);
    mpreal rootbig   = one / rootsfmin;
    printf("  rootsfmin   = "); _printnum_full(rootsfmin); printf("\n");
    printf("  rootbig     = "); _printnum_full(rootbig);   printf("\n");
    printf("\n");

    // ---- 2. Test with dense matrix ----
    int test_sizes[] = {3, 5, 10};
    int num_tests = 3;

    for (int ti = 0; ti < num_tests; ti++) {
        mplapackint n = test_sizes[ti];
        mplapackint m = n;
        mplapackint lda = m;
        mplapackint ldv = n;

        printf("============================================\n");
        printf("  Test M=%d, N=%d (type 4, DENSE)\n", (int)m, (int)n);
        printf("============================================\n");

        mpcomplex *a     = new mpcomplex[lda * n];
        mpcomplex *a_sav = new mpcomplex[lda * n];
        mpreal    *sva   = new mpreal[n];
        mpreal    *sv_target = new mpreal[n];
        mpcomplex *v     = new mpcomplex[ldv * n];
        mplapackint lwork  = m + n + 10;
        mplapackint lrwork = max((mplapackint)6, m + n);
        mpcomplex *cwork = new mpcomplex[lwork];
        mpreal    *rwork = new mpreal[lrwork];

        // Step 1: Create diagonal matrix with type 4 singular values
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                a[i + j * lda] = mpcomplex(zero, zero);
            }
            // Evenly spaced near underflow: sv = sfmin * (n - j) / n
            // (this gives values between sfmin/n and sfmin)
            // Actually, zdrvbd uses: sv(i) = ulp + (1 - ulp) * (n-i)/(n-1) * sfmin
            // We use a simpler version:
            sv_target[j] = sfmin * mpreal(n - j);
            a[j + j * lda] = mpcomplex(sv_target[j], zero);
        }

        printf("  Target singular values (before transforms):\n");
        for (int i = 0; i < n; i++) {
            printf("    sv[%d] = ", i); _printnum_full(sv_target[i]); printf("\n");
        }

        // Step 2: Apply multiple Householder reflections to make it dense
        // A <- H_L3 * H_L2 * H_L1 * diag(sv) * H_R1 * H_R2 * H_R3
        int num_transforms = 5;
        for (int t = 0; t < num_transforms; t++) {
            apply_random_left_householder(m, n, a, lda, 12345 + t * 1000 + ti * 100);
            apply_random_right_householder(m, n, a, lda, 67890 + t * 1000 + ti * 100);
        }

        // Save A
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                a_sav[i + j * lda] = a[i + j * lda];

        // Step 3: Compute column norms using safe method
        printf("  Column norms of dense A (safe computation):\n");
        for (int j = 0; j < n; j++) {
            // Use Rnrm2-like safe computation
            mpreal scale_c = zero, ssq_c = one;
            for (int i = 0; i < m; i++) {
                mpreal re = abs(a[i + j * lda].real());
                mpreal im = abs(a[i + j * lda].imag());
                if (re != zero) {
                    if (scale_c < re) {
                        ssq_c = one + ssq_c * (scale_c / re) * (scale_c / re);
                        scale_c = re;
                    } else {
                        ssq_c += (re / scale_c) * (re / scale_c);
                    }
                }
                if (im != zero) {
                    if (scale_c < im) {
                        ssq_c = one + ssq_c * (scale_c / im) * (scale_c / im);
                        scale_c = im;
                    } else {
                        ssq_c += (im / scale_c) * (im / scale_c);
                    }
                }
            }
            mpreal cnorm = scale_c * sqrt(ssq_c);
            printf("    ||A(:,%d)|| = ", j); _printnum(cnorm);
            printf("  (< rootsfmin? %s)\n",
                   (cnorm < rootsfmin) ? "YES" : "no");
        }

        mpreal anorm_f = safe_frobenius_norm_complex(m, n, a, lda);
        printf("  ||A||_F (safe) = "); _printnum_full(anorm_f); printf("\n");

        // Step 4: Call Cgesvj
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                v[i + j * ldv] = (i == j) ? mpcomplex(one, zero)
                                           : mpcomplex(zero, zero);

        mplapackint info = 0;
        printf("\n  Calling Cgesvj('G', 'U', 'V', %d, %d, ...)...\n",
               (int)m, (int)n);
        Cgesvj("G", "U", "V", m, n, a, lda, sva, 0, v, ldv,
               cwork, lwork, rwork, lrwork, info);

        printf("  INFO    = %d\n", (int)info);
        printf("  SCALE   = RWORK(1)  = "); _printnum_full(rwork[0]); printf("\n");
        printf("  NRANK   = RWORK(2)  = "); _printnum(rwork[1]); printf("\n");
        printf("  NRANK>S = RWORK(3)  = "); _printnum(rwork[2]); printf("\n");
        printf("  NSWEEP  = RWORK(4)  = "); _printnum(rwork[3]); printf("\n");
        printf("  maxcos  = RWORK(5)  = "); _printnum(rwork[4]); printf("\n");
        printf("  maxsin  = RWORK(6)  = "); _printnum(rwork[5]); printf("\n");

        // Also check CWORK(1) in case SCALE is stored there
        printf("  CWORK(1) = ("); _printnum(cwork[0].real());
        printf(", "); _printnum(cwork[0].imag()); printf(")\n");

        mpreal scale = rwork[0];

        // Step 5: Computed singular values
        printf("\n  Computed singular values:\n");
        for (int i = 0; i < n; i++) {
            mpreal sv_comp = scale * sva[i];
            mpreal sv_tgt  = sv_target[i];  // sorted descending
            mpreal relerr;
            if (sv_tgt != zero) {
                relerr = abs(sv_comp - sv_tgt) / sv_tgt;
            } else {
                relerr = abs(sv_comp);
            }
            printf("    sva[%d] = ", i); _printnum(sva[i]);
            printf("  scaled = "); _printnum(sv_comp);
            printf("  target = "); _printnum(sv_tgt);
            printf("  relerr = "); _printnum(relerr);
            printf("\n");
        }

        // Step 6: Compute residual ||A - U diag(S) V^H|| using safe norm
        // Build R = A_sav - U * diag(scale*sva) * V^H
        mpcomplex *resid = new mpcomplex[lda * n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                resid[i + j * lda] = a_sav[i + j * lda];
                for (int k = 0; k < n; k++) {
                    mpcomplex u_ik = a[i + k * lda];
                    mpcomplex v_jk = v[j + k * ldv];
                    mpreal    s_k  = scale * sva[k];
                    resid[i + j * lda] -= u_ik * mpcomplex(s_k, zero) * conj(v_jk);
                }
            }
        }

        mpreal resid_f = safe_frobenius_norm_complex(m, n, resid, lda);
        mpreal ulp = eps;
        mpreal mn_max = mpreal(max(m, n));
        mpreal ratio;
        if (anorm_f > zero) {
            ratio = resid_f / (anorm_f * mn_max * ulp);
        } else {
            if (resid_f > zero) {
                ratio = resid_f / (mn_max * ulp); // avoid 0/0
                printf("  WARNING: ||A||=0 but ||resid||!=0\n");
            } else {
                ratio = zero;
            }
        }

        printf("\n  --- Test 15 (decomposition accuracy) ---\n");
        printf("  ||A||_F (safe)        = "); _printnum_full(anorm_f); printf("\n");
        printf("  ||A - USVT||_F (safe) = "); _printnum_full(resid_f); printf("\n");
        printf("  ratio = resid / (||A|| * N * ulp) = "); _printnum_full(ratio); printf("\n");
        printf("  threshold = 50.0\n");
        printf("  PASS? %s\n", (ratio < 50.0) ? "YES" : "*** FAIL ***");

        // Step 7: Check 1/(2*eps)
        mpreal inv_2eps = one / (2.0 * eps);
        printf("  1/(2*eps) = "); _printnum_full(inv_2eps); printf("\n");
        mpreal ratio_diff = abs(ratio - inv_2eps);
        printf("  |ratio - 1/(2*eps)| = "); _printnum(ratio_diff); printf("\n");
        if (ratio_diff < one) {
            printf("  *** ratio  1/(2*eps) : SYSTEMATIC SCALING FAILURE ***\n");
        }

        // Step 8: U orthogonality
        mpreal orth_u_f = zero;
        {
            mpreal orth_scale = zero, orth_ssq = one;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    mpcomplex dot = mpcomplex(zero, zero);
                    for (int k = 0; k < m; k++) {
                        dot += conj(a[k + i * lda]) * a[k + j * lda];
                    }
                    mpcomplex diff = dot - ((i == j) ? mpcomplex(one, zero)
                                                     : mpcomplex(zero, zero));
                    mpreal re = abs(diff.real());
                    mpreal im = abs(diff.imag());
                    if (re != zero) {
                        if (orth_scale < re) {
                            orth_ssq = one + orth_ssq * (orth_scale/re)*(orth_scale/re);
                            orth_scale = re;
                        } else {
                            orth_ssq += (re/orth_scale)*(re/orth_scale);
                        }
                    }
                    if (im != zero) {
                        if (orth_scale < im) {
                            orth_ssq = one + orth_ssq * (orth_scale/im)*(orth_scale/im);
                            orth_scale = im;
                        } else {
                            orth_ssq += (im/orth_scale)*(im/orth_scale);
                        }
                    }
                }
            }
            orth_u_f = orth_scale * sqrt(orth_ssq);
        }
        mpreal ratio_u = orth_u_f / (mpreal(m) * ulp);
        printf("\n  --- Test 16 (U orthogonality) ---\n");
        printf("  ||I - U^H U||_F (safe) = "); _printnum_full(orth_u_f); printf("\n");
        printf("  ratio                  = "); _printnum_full(ratio_u); printf("\n");
        printf("  PASS? %s\n", (ratio_u < 50.0) ? "YES" : "*** FAIL ***");

        // Step 9: What if we ignore SCALE?
        printf("\n  --- What if SCALE is ignored? ---\n");
        if (scale != one) {
            mpcomplex *resid2 = new mpcomplex[lda * n];
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < m; i++) {
                    resid2[i + j * lda] = a_sav[i + j * lda];
                    for (int k = 0; k < n; k++) {
                        mpcomplex u_ik = a[i + k * lda];
                        mpcomplex v_jk = v[j + k * ldv];
                        mpreal    s_k  = sva[k]; // NO scale!
                        resid2[i + j * lda] -= u_ik * mpcomplex(s_k, zero) * conj(v_jk);
                    }
                }
            }
            mpreal resid2_f = safe_frobenius_norm_complex(m, n, resid2, lda);
            mpreal ratio2;
            if (anorm_f > zero) {
                ratio2 = resid2_f / (anorm_f * mn_max * ulp);
            } else {
                ratio2 = zero;
            }
            printf("  SCALE = "); _printnum(scale); printf(" (not 1.0!)\n");
            printf("  ratio WITHOUT scale = "); _printnum_full(ratio2); printf("\n");
            printf("  ratio WITH scale    = "); _printnum_full(ratio);  printf("\n");
            delete[] resid2;
        } else {
            printf("  SCALE = 1.0 (no difference)\n");
        }

        delete[] resid;
        delete[] a;
        delete[] a_sav;
        delete[] sva;
        delete[] sv_target;
        delete[] v;
        delete[] cwork;
        delete[] rwork;

        printf("\n");
    }

    return 0;
}
