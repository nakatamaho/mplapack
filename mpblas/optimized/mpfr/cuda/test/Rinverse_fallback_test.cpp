// Test the public inverse dispatchers in libmplapack_mpfr_opt_cuda, with
// and without the fixed-precision MPFR CPU hooks. No CUDA device is needed.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <sys/wait.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_mpfr_cuda.h>
#include "../../openmp/Rfreal_state.h"

namespace {
using I = mplapackint;
enum Routine { trtri, lauum, getri, potri };
constexpr const char *names[] = {"trtri", "lauum", "getri", "potri"};

template <class T> struct Api;
#define INVERSE_API(T, P) \
    template <> struct Api<T> { \
        static constexpr char prefix = #P[0]; \
        template <class F> static void public_call(F f) { f(P##trtri, P##lauum, P##getri, P##potri); } \
        template <class F> static bool cuda(F f) { return f(P##trtri_mpfr_cuda, P##lauum_mpfr_cuda, P##getri_mpfr_cuda, P##potri_mpfr_cuda); } \
    }
INVERSE_API(mpfr_class, R);
INVERSE_API(mpc_class, C);
#undef INVERSE_API

void set_env(const char *key, const char *value)
{
    if (setenv(key, value, 1) != 0) {
        std::perror(key);
        std::exit(1);
    }
}

void clear_env()
{
    // Collect names first: unsetenv can move entries in environ.
    std::vector<std::string> keys;
    for (char **entry = environ; *entry; ++entry) {
        const std::string var(*entry);
        if (var.compare(0, 9, "MPLAPACK_") == 0 || var.compare(0, 7, "MPFRXX_") == 0 || var.compare(0, 6, "MPCXX_") == 0)
            keys.push_back(var.substr(0, var.find('=')));
    }
    for (const auto &key : keys) {
        if (unsetenv(key.c_str()) != 0) {
            std::perror(key.c_str());
            std::exit(1);
        }
    }
}

mpfr_class scalar(double re, double, mpfr_class *) { return mpfr_class(re); }
mpc_class scalar(double re, double im, mpc_class *) { return mpc_class(re, im); }

// Store the full significand, plus precision, in the freal-off process.
// Check exact equality (including the sign of zero), not a residual or norm.
bool match(FILE *ref, bool record, mpfr_srcptr value, mpfr_ptr expected)
{
    const long precision = mpfr_get_prec(value);
    if (record)
        return std::fprintf(ref, "%ld ", precision) > 0 && mpfr_out_str(ref, 16, 0, value, MPFR_RNDN) != 0 && std::fputc('\n', ref) != EOF;
    long ref_precision;
    mpfr_set_prec(expected, precision);
    return std::fscanf(ref, "%ld", &ref_precision) == 1 && ref_precision == precision &&
           mpfr_inp_str(expected, ref, 16, MPFR_RNDN) != 0 && mpfr_equal_p(value, expected) &&
           mpfr_signbit(value) == mpfr_signbit(expected);
}

bool match(FILE *ref, bool record, const mpfr_class &value, mpfr_ptr expected)
{
    return match(ref, record, value.get_mpfr_t(), expected);
}

bool match(FILE *ref, bool record, const mpc_class &value, mpfr_ptr expected)
{
    return match(ref, record, mpc_realref(value.mpc_data()), expected) && match(ref, record, mpc_imagref(value.mpc_data()), expected);
}

bool identical(mpfr_srcptr x, mpfr_srcptr y)
{
    return mpfr_get_prec(x) == mpfr_get_prec(y) && mpfr_equal_p(x, y) && mpfr_signbit(x) == mpfr_signbit(y);
}

bool identical(const mpfr_class &x, const mpfr_class &y)
{
    return identical(x.get_mpfr_t(), y.get_mpfr_t());
}

bool identical(const mpc_class &x, const mpc_class &y)
{
    return identical(mpc_realref(x.mpc_data()), mpc_realref(y.mpc_data())) && identical(mpc_imagref(x.mpc_data()), mpc_imagref(y.mpc_data()));
}

mpfr_class conjugate(const mpfr_class &x) { return x; }
mpc_class conjugate(const mpc_class &x) { return conj(x); }

// Reconstruct the original matrix and expand the stored inverse without
// using any LAPACK or BLAS routines. Use working-precision dot products.
template <class T>
bool inverse_residual(Routine routine, bool upper, bool unit, long precision, I n, I lda, const std::vector<T> &input,
                      const std::vector<I> &piv, const std::vector<T> &output)
{
    const T zero(0.0), one(1.0);
    std::vector<T> orig(n * n, zero), inv(n * n, zero);
    for (I j = 0; j < n; ++j) {
        for (I i = 0; i < n; ++i) {
            const I ij = i + j * n;
            const bool stored = upper ? i <= j : i >= j;
            if (routine == getri) {
                // Packed unit lower L and upper U.
                for (I k = 0; k < std::min(i, j); ++k)
                    orig[ij] += input[i + k * lda] * input[k + j * lda];
                if (i <= j)
                    orig[ij] += input[i + j * lda];
                else
                    orig[ij] += input[i + j * lda] * input[j + j * lda];
                inv[ij] = output[i + j * lda];
            } else if (routine == potri) {
                for (I k = 0; k <= std::min(i, j); ++k) {
                    if (upper)
                        orig[ij] += conjugate(input[k + i * lda]) * input[k + j * lda];
                    else
                        orig[ij] += input[i + k * lda] * conjugate(input[j + k * lda]);
                }
                inv[ij] = stored ? output[i + j * lda] : conjugate(output[j + i * lda]);
            } else if (stored) {
                orig[ij] = unit && i == j ? one : input[i + j * lda];
                inv[ij] = unit && i == j ? one : output[i + j * lda];
            }
        }
    }
    if (routine == getri) {
        // Undo GETRF's row swaps in reverse order: A = P*L*U.
        for (I i = n; i-- > 0;)
            if (piv[i] != i + 1)
                for (I j = 0; j < n; ++j)
                    std::swap(orig[i + j * n], orig[piv[i] - 1 + j * n]);
    }

    mpfr_class bound;
    mpfr_set_ui_2exp(bound.get_mpfr_t(), 16 * n, -precision, MPFR_RNDN);
    for (I j = 0; j < n; ++j) {
        for (I i = 0; i < n; ++i) {
            T error = i == j ? T(-1.0) : zero;
            for (I k = 0; k < n; ++k)
                error += orig[i + k * n] * inv[k + j * n];
            const mpfr_class magnitude = abs(error);
            if (!mpfr_lessequal_p(magnitude.get_mpfr_t(), bound.get_mpfr_t()))
                return false;
        }
    }
    return true;
}

template <class T>
bool run_case(FILE *ref, bool record, bool freal_available, long precision, I n, Routine routine, const char *uplo, const char *diag)
{
    using A = Api<T>;
    const I lda = n + 3;
    const char getri_name[] = {A::prefix, 'g', 'e', 't', 'r', 'i', '\0'};
    const I lwork = routine == getri ? std::max<I>(1, n * iMlaenv(1, getri_name, " ", n, -1, -1, -1)) : 0;
    const auto make = [](double re, double im) { return scalar(re, im, static_cast<T *>(nullptr)); };
    const T sentinel = make(-99, 33);
    std::vector<T> a(lda * std::max<I>(1, n), sentinel), work(lwork, sentinel);
    std::vector<I> piv(n);
    for (I j = 0; j < n; ++j) {
        // GETRI takes packed L/U factors, not the unfactored matrix.
        piv[j] = j == 0 && n > 1 ? 2 : j + 1;
        for (I i = 0; i < n; ++i) {
            if (routine != getri && (*uplo == 'U' ? i > j : i < j))
                continue;
            if (i == j && *diag == 'U')
                continue; // A unit-diagonal inverse must not touch the sentinel.
            const double re = i == j ? 3 + j % 5 : (1 + (i + 3 * j) % 7) / 64.0;
            const double im = i == j && routine == potri ? 0 : (1 + (3 * i + j) % 5) / 64.0;
            a[i + j * lda] = make(re, im);
        }
    }
    I info = -123;
    const auto invoke = [&](auto tri, auto lau, auto get, auto pot) {
        switch (routine) {
        case trtri: return tri(uplo, diag, n, a.data(), lda, info);
        case lauum: return lau(uplo, n, a.data(), lda, info);
        case getri: return get(n, a.data(), lda, piv.data(), work.data(), lwork, info);
        default: return pot(uplo, n, a.data(), lda, info);
        }
    };
    const auto fail = [&](const char *what, long detail) {
        std::fprintf(stderr, "%c%s p=%ld n=%ld uplo=%s diag=%s freal=%s: %s (%ld)\n",
                     A::prefix, names[routine], precision, long(n), uplo, diag, record ? "off" : "on", what, detail);
        return false;
    };

    const auto a_input = a, work_input = work;
    const auto piv_input = piv;
    // MPLAPACK_MPFR_CUDA=0 disables all CUDA entries (including BLAS).
    // The inverse probe declines before taking the GPU lock or setting up a device.
    if (!record) {
        const I info_input = info;
        if (A::cuda(invoke) || info != info_input)
            return fail("CUDA entry did not decline cleanly", info);
        for (size_t i = 0; i < a.size(); ++i)
            if (!identical(a[i], a_input[i]))
                return fail("declined CUDA entry changed matrix at index", i);
        for (size_t i = 0; i < work.size(); ++i)
            if (!identical(work[i], work_input[i]))
                return fail("declined CUDA entry changed workspace at index", i);
        if (piv != piv_input)
            return fail("declined CUDA entry changed pivots", 0);
    }
    const auto before = mplapack_freal_done.load();
    A::public_call(invoke); // Public LAPACK symbol, supplied only by the CUDA library.
    const auto done = mplapack_freal_done.load() - before;
    if (info != 0)
        return fail("unexpected info", info);
    const unsigned long hooks = !record && freal_available && precision <= 512 && n >= 7 ? (routine == potri ? 2 : 1) : 0;
    // Only n=7 isolates the inverse hooks: the minimum work (64) admits
    // them but excludes individual BLAS calls. POTRI must run both TRTRI
    // and LAUUM. At n=67, blocked GETRI can meet the count via BLAS alone.
    if (done < hooks || ((routine != getri || n <= 7 || hooks == 0) && done != hooks))
        return fail("unexpected number of freal completions", done);

    const bool upper = *uplo == 'U', unit = *diag == 'U';
    for (I j = 0; j < std::max<I>(1, n); ++j) {
        for (I i = 0; i < lda; ++i) {
            const bool untouched = i >= n || (routine != getri && (upper ? i > j : i < j)) || (routine == trtri && unit && i == j);
            if (untouched && !identical(a[i + j * lda], sentinel))
                return fail("public routine changed sentinel at index", i + j * lda);
        }
    }
    if (!record && precision <= 512 && n >= 7 && routine != lauum &&
        !inverse_residual(routine, upper, unit, precision, n, lda, a_input, piv_input, a))
        return fail("inverse residual exceeds 16*n*eps", 0);

    mpfr_class expected;
    for (size_t i = 0; i < a.size(); ++i)
        if (!match(ref, record, a[i], expected.get_mpfr_t()))
            return fail("matrix differs from freal-off result at index", i);
    for (size_t i = 0; i < work.size(); ++i)
        if (!match(ref, record, work[i], expected.get_mpfr_t()))
            return fail("workspace differs from freal-off result at index", i);
    return true;
}

bool run_suite(FILE *ref, bool record, bool freal_available, int &checks)
{
    set_env("MPLAPACK_MPFR_CPU_FREAL", record ? "0" : "1");
    for (long precision : {256L, 512L, 4096L}) {
        mpfrxx::set_default_precision_bits(precision);
        mpfrxx::set_default_mpc_precision_bits(precision);
        // Pin defaults after wrapper initialization in each process.
        mpfrxx::set_default_rounding_mode(MPFR_RNDN);
        mpfrxx::set_default_mpc_rounding_mode(MPFR_RNDN);
        if (mpfr_set_emin(MPFR_EMIN_DEFAULT) != 0 || mpfr_set_emax(MPFR_EMAX_DEFAULT) != 0) {
            std::fputs("could not reset MPFR exponent range\n", stderr);
            return false;
        }
        // 67 crosses the 64-element TRTRI, LAUUM and GETRI block sizes.
        for (I n : {0, 1, 7, 67}) {
            if (precision == 4096 && n == 67)
                continue;
            const auto test = [&](Routine routine, const char *uplo, const char *diag) {
                checks += 2;
                return run_case<mpfr_class>(ref, record, freal_available, precision, n, routine, uplo, diag) &&
                       run_case<mpc_class>(ref, record, freal_available, precision, n, routine, uplo, diag);
            };
            if (!test(getri, "U", "N"))
                return false;
            for (const char *uplo : {"U", "L"})
                if (!test(trtri, uplo, "N") || !test(trtri, uplo, "U") || !test(lauum, uplo, "N") || !test(potri, uplo, "N"))
                    return false;
        }
    }
    return true;
}
} // namespace

int main()
{
    clear_env();
    set_env("MPLAPACK_MPFR_CUDA", "0");
    set_env("MPLAPACK_MPFR_LAPACK_CUDA", "1");
    set_env("MPLAPACK_MPFR_LAPACK_INVERSE_CUDA_MIN", "0");
    set_env("MPLAPACK_MPFR_CPU_FREAL_MIN", "64");
    set_env("MPLAPACK_MPFR_CPU_FREAL_L12_PREC", "512");
    set_env("MPLAPACK_MPFR_GEMM_WINOGRAD_CUTOFF", "0");
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    bool freal_available = true;
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    freal_available = __builtin_cpu_supports("bmi2") && __builtin_cpu_supports("adx");
#endif
    if (!freal_available)
        std::puts("CPU lacks BMI2/ADX; checking the freal-decline paths only");

    // The freal enable flag is cached at first use. Fork before either side
    // calls LAPACK or BLAS, so both runs use the same public fallback path but
    // have independent caches. The child writes the reference to a temp file.
    FILE *ref = std::tmpfile();
    if (!ref) {
        std::perror("tmpfile");
        return 1;
    }
    const pid_t pid = fork();
    int checks = 0;
    if (pid == 0) {
        const bool ok = run_suite(ref, true, freal_available, checks);
        _exit(std::fclose(ref) == 0 && ok ? 0 : 1);
    }
    if (pid < 0) {
        std::perror("fork");
        return 1;
    }
    int status;
    if (waitpid(pid, &status, 0) != pid || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        std::fputs("freal-off reference process failed\n", stderr);
        return 1;
    }
    if (std::fseek(ref, 0, SEEK_SET) != 0) {
        std::perror("fseek");
        return 1;
    }
    const bool ok = run_suite(ref, false, freal_available, checks);
    std::fclose(ref);
    if (!ok)
        return 1;
    std::printf("inverse CUDA CPU fallback: %d cases matched bit for bit\n", checks);
    return 0;
}
