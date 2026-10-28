// Factorization regression gates. No library routines are replaced by this test.
#include <mpblas.h>
#include <mplapack.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <set>
#include <type_traits>
#include <vector>
#ifdef MPLAPACK_FACTOR_CUDA
#include <cuda_runtime_api.h>
#endif

void Rgetrf_reference_cpu(INTEGER, INTEGER, REAL *, INTEGER, INTEGER *, INTEGER &);
void Cgetrf_reference_cpu(INTEGER, INTEGER, COMPLEX *, INTEGER, INTEGER *, INTEGER &);
void Rgetrf2_reference_cpu(INTEGER, INTEGER, REAL *, INTEGER, INTEGER *, INTEGER &);
void Cgetrf2_reference_cpu(INTEGER, INTEGER, COMPLEX *, INTEGER, INTEGER *, INTEGER &);
void Rpotrf_reference_cpu(const char *, INTEGER, REAL *, INTEGER, INTEGER &);
void Cpotrf_reference_cpu(const char *, INTEGER, COMPLEX *, INTEGER, INTEGER &);

namespace {
int checks = 0, failures = 0;
const char *routine = "";
char triangle = '-';
bool fast;
REAL worst = 0;
std::mt19937 rng(20261101);
REAL random_real() {
    REAL x = 0;
    // Exercise all limbs, not just the leading double.
    for (int k = 0; k < 10; ++k)
        x += REAL(std::ldexp(double(int(rng() % 65536) - 32768), -16 - 48*k));
    return x;
}
REAL magnitude(const REAL &x) { return abs(x); }
REAL magnitude(const COMPLEX &x) { return abs(x.real()) + abs(x.imag()); }
REAL adjoint(const REAL &x) { return x; }
COMPLEX adjoint(const COMPLEX &x) { return conj(x); }
REAL random_scalar(REAL *) { return random_real(); }
COMPLEX random_scalar(COMPLEX *) { return COMPLEX(random_real(), random_real()); }
bool same(const REAL &a, const REAL &b) {
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    // Object bytes contain allocation pointers. Compare the exact numeric
    // representation instead (including signed zero and precision).
    return mpfr_get_prec(a.get_mpfr_t()) == mpfr_get_prec(b.get_mpfr_t()) &&
        mpfr_equal_p(a.get_mpfr_t(), b.get_mpfr_t()) &&
        mpfr_signbit(a.get_mpfr_t()) == mpfr_signbit(b.get_mpfr_t());
#elif defined(MPLAPACK_BUILD_WITH_GMP)
    return mpf_cmp(a.get_mpf_t(), b.get_mpf_t()) == 0;
#elif defined(MPLAPACK_BUILD_WITH_BINARY80)
    // x87 extended has six indeterminate padding bytes.
    return a == b && std::signbit(a) == std::signbit(b);
#else
    return std::memcmp(&a, &b, sizeof a) == 0;
#endif
}
bool same(const COMPLEX &a, const COMPLEX &b) {
    return same(REAL(a.real()), REAL(b.real())) && same(REAL(a.imag()), REAL(b.imag()));
}
void product(INTEGER m, INTEGER n, INTEGER k, const REAL *l, const REAL *u, REAL *p) {
    Rgemm("N", "N", m, n, k, REAL(1), const_cast<REAL *>(l), m,
          const_cast<REAL *>(u), k, REAL(0), p, m);
}
void product(INTEGER m, INTEGER n, INTEGER k, const COMPLEX *l, const COMPLEX *u, COMPLEX *p) {
    Cgemm("N", "N", m, n, k, COMPLEX(1), const_cast<COMPLEX *>(l), m,
          const_cast<COMPLEX *>(u), k, COMPLEX(0), p, m);
}
void lu(INTEGER m, INTEGER n, REAL *a, INTEGER ld, INTEGER *p, INTEGER &i, bool ref) {
    if (ref) Rgetrf_reference_cpu(m,n,a,ld,p,i); else Rgetrf(m,n,a,ld,p,i);
}
void lu(INTEGER m, INTEGER n, COMPLEX *a, INTEGER ld, INTEGER *p, INTEGER &i, bool ref) {
    if (ref) Cgetrf_reference_cpu(m,n,a,ld,p,i); else Cgetrf(m,n,a,ld,p,i);
}
void panel(INTEGER m, INTEGER n, REAL *a, INTEGER ld, INTEGER *p, INTEGER &i, bool ref) {
    if (ref) Rgetrf2_reference_cpu(m,n,a,ld,p,i); else Rgetrf2(m,n,a,ld,p,i);
}
void panel(INTEGER m, INTEGER n, COMPLEX *a, INTEGER ld, INTEGER *p, INTEGER &i, bool ref) {
    if (ref) Cgetrf2_reference_cpu(m,n,a,ld,p,i); else Cgetrf2(m,n,a,ld,p,i);
}
void chol(const char *u, INTEGER n, REAL *a, INTEGER ld, INTEGER &i, bool ref) {
    if (ref) Rpotrf_reference_cpu(u,n,a,ld,i); else Rpotrf(u,n,a,ld,i);
}
void chol(const char *u, INTEGER n, COMPLEX *a, INTEGER ld, INTEGER &i, bool ref) {
    if (ref) Cpotrf_reference_cpu(u,n,a,ld,i); else Cpotrf(u,n,a,ld,i);
}
void check(bool ok, const char *what, INTEGER m, INTEGER n, int kind) {
    ++checks;
    if (!ok) { ++failures; std::printf("FAIL %s(%c) %s %ldx%ld kind=%d\n", routine, triangle, what, long(m), long(n), kind); }
}
template<class T> void identical(const std::vector<T> &r, const std::vector<T> &g,
                                INTEGER m, INTEGER n, int kind) {
    bool ok = true;
    for (size_t t=0; t<r.size(); ++t) ok = same(r[t],g[t]) && ok;
    check(ok,"bit identity",m,n,kind);
}
template<class T> void residual(const std::vector<T> &a, const std::vector<T> &p,
                               INTEGER m, INTEGER n, INTEGER ld, int kind) {
    REAL anorm=0, error=0;
    for (INTEGER j=0;j<n;++j) {
        REAL ac=0, ec=0;
        for (INTEGER i=0;i<m;++i) { ac+=magnitude(a[i+j*ld]); ec+=magnitude(a[i+j*ld]-p[i+j*m]); }
        if (ac>anorm) anorm=ac;
        if (ec>error) error=ec;
    }
    REAL ratio = anorm == REAL(0) ? error : REAL(error/(REAL(n)*anorm*Rlamch("E")));
    if (ratio>worst) worst=ratio;
    check(ratio<=REAL(30),"backward error",m,n,kind);
}
template<class T> void test_lu(INTEGER m, INTEGER n, int kind, bool direct_panel = false) {
    routine = std::is_same<T, REAL>::value ? "Rgetrf" : "Cgetrf";
    triangle = '-';
    INTEGER ld=m+3, k=std::min(m,n), ir=-99, ig=-98;
    std::vector<T> a(ld*n,T(19));
    for (INTEGER j=0;j<n;++j) for (INTEGER i=0;i<m;++i) a[i+j*ld]=random_scalar((T*)nullptr);
    // zero column; exactly dependent columns; forced swaps at every step;
    // tied maxima (first-index semantics).
    if (kind==1) for (INTEGER i=0;i<m;++i) a[i+(k/2)*ld]=T(0);
    if (kind==2) for (INTEGER i=0;i<m;++i) a[i+(n-1)*ld]=a[i];
    if (kind>=3) {
        for (INTEGER j=0;j<n;++j) for (INTEGER i=0;i<m;++i) a[i+j*ld]=T(0);
        for (INTEGER j=0;j<k;++j) a[m-1-j+j*ld]=T(2);
        if (kind==4) for (INTEGER i=0;i<m;++i) a[i]=T(2);
    }
    auto r=a, g=a;
    std::vector<INTEGER> pr(k,-1), pg(k,-2);
    // The oracle must always use the blocked reference path.
    if (direct_panel) {
        panel(m,n,r.data(),ld,pr.data(),ir,true);
        panel(m,n,g.data(),ld,pg.data(),ig,false);
    } else {
        lu(m,n,r.data(),ld,pr.data(),ir,true);
        lu(m,n,g.data(),ld,pg.data(),ig,false);
    }
    check(ir==ig,"LU info",m,n,kind);
    if (!fast || kind!=4) check(std::memcmp(pr.data(),pg.data(),k*sizeof(INTEGER))==0,"pivots",m,n,kind);
    if (!fast) identical(r,g,m,n,kind);
    std::vector<T> l(m*k,T(0)), u(k*n,T(0)), p(m*n,T(0));
    for (INTEGER j=0;j<k;++j) for (INTEGER i=j;i<m;++i) l[i+j*m]=i==j?T(1):g[i+j*ld];
    for (INTEGER j=0;j<n;++j) for (INTEGER i=0;i<std::min(k,j+1);++i) u[i+j*k]=g[i+j*ld];
    for (INTEGER i=0;i<k;++i) {
        if (pg[i]<1 || pg[i]>m) { check(false,"pivot bounds",m,n,kind); return; }
        for (INTEGER j=0;j<n;++j) std::swap(a[i+j*ld],a[pg[i]-1+j*ld]);
    }
    product(m,n,k,l.data(),u.data(),p.data());
    residual(a,p,m,n,ld,kind);
}
// Complex -one has a negative-zero imaginary component, unlike COMPLEX(-1).
// Exercise both the new recursive panel and the serial blocked dispatcher.
void test_signed_zeros() {
    const char *old = std::getenv("MPLAPACK_LAPACK_LOOKAHEAD");
    const bool had_old = old != nullptr;
    const std::string saved = old ? old : "";
    setenv("MPLAPACK_LAPACK_LOOKAHEAD", "0", 1);
    for (bool direct : {true, false}) {
        if (!direct && fast) continue; // FAST getrf may round differently.
        const INTEGER m = direct ? 2 : 65, n = 65, ld = m + 3;
        routine = direct ? "Cgetrf2" : "Cgetrf";
        triangle = '-';
        for (int pa=0; pa<4; ++pa) for (int la=0; la<4; ++la)
        for (int ba=0; ba<4; ++ba) for (int ca=0; ca<4; ++ca) {
            auto zero = [](int flags) {
                return COMPLEX(REAL(flags&1 ? -0.0 : 0.0), REAL(flags&2 ? -0.0 : 0.0));
            };
            const COMPLEX pivot(REAL(pa&1 ? -1.0 : 1.0), REAL(pa&2 ? -0.0 : 0.0));
            const COMPLEX upper(REAL(ba&1 ? -1.0 : 1.0), REAL(ba&2 ? -0.0 : 0.0));
            std::vector<COMPLEX> a(ld*n, COMPLEX(19));
            for (INTEGER j=0; j<n; ++j) for (INTEGER i=0; i<m; ++i) a[i+j*ld]=COMPLEX(0);
            if (direct) {
                a[0]=pivot; a[1]=zero(la);
                for (INTEGER j=1; j<n; ++j) { a[j*ld]=upper; a[1+j*ld]=zero(ca); }
            } else {
                for (INTEGER j=0; j<m-1; ++j) a[j+j*ld]=COMPLEX(1);
                a[m-2+(m-2)*ld]=pivot; a[m-1+(m-2)*ld]=zero(la);
                a[m-2+(n-1)*ld]=upper; a[m-1+(n-1)*ld]=zero(ca);
            }
            auto r=a, g=a;
            std::vector<INTEGER> pr(m), pg(m);
            INTEGER ir=0, ig=0;
            if (direct) {
                panel(m,n,r.data(),ld,pr.data(),ir,true);
                panel(m,n,g.data(),ld,pg.data(),ig,false);
            } else {
                lu(m,n,r.data(),ld,pr.data(),ir,true);
                lu(m,n,g.data(),ld,pg.data(),ig,false);
            }
            const int kind=pa+4*la+16*ba+64*ca;
            check(ir==ig,"signed-zero info",m,n,kind);
            check(std::memcmp(pr.data(),pg.data(),m*sizeof(INTEGER))==0,"signed-zero pivots",m,n,kind);
            identical(r,g,m,n,kind);
        }
    }
    if (had_old) setenv("MPLAPACK_LAPACK_LOOKAHEAD", saved.c_str(), 1);
    else unsetenv("MPLAPACK_LAPACK_LOOKAHEAD");
}
template<class T> void test_chol(INTEGER n, const char *up, int bad) {
    routine = std::is_same<T, REAL>::value ? "Rpotrf" : "Cpotrf";
    triangle = *up;
    INTEGER ld=n+3, ir=-99, ig=-98;
    std::vector<T> a(ld*n,T(19));
    for (INTEGER j=0;j<n;++j) for (INTEGER i=0;i<=j;++i) {
        T x = i==j ? T(2*n) : random_scalar((T*)nullptr);
        a[i+j*ld]=x; a[j+i*ld]=adjoint(x);
    }
    // Isolate the bad pivot so the expected info is known exactly.
    INTEGER fail_at=n/2;
    if (bad) {
        for (INTEGER j=0;j<n;++j) a[fail_at+j*ld]=a[j+fail_at*ld]=T(0);
        a[fail_at+fail_at*ld]=T(-1);
    }
    auto r=a,g=a;
    chol(up,n,r.data(),ld,ir,true); chol(up,n,g.data(),ld,ig,false);
    check(ir==ig && ig==(bad?fail_at+1:0),"Cholesky info",n,n,bad);
    if (!fast) identical(r,g,n,n,bad);
    if (bad) return;
    std::vector<T> l(n*n,T(0)),u(n*n,T(0)),p(n*n,T(0));
    for (INTEGER j=0;j<n;++j) for (INTEGER i=j;i<n;++i) {
        T x = *up=='L'?g[i+j*ld]:adjoint(g[j+i*ld]);
        l[i+j*n]=x; u[j+i*n]=adjoint(x);
    }
    product(n,n,n,l.data(),u.data(),p.data()); residual(a,p,n,n,ld,bad);
}
}
int main() {
#ifndef MPLAPACK_FACTOR_CUDA
    // CPU: exercise tasks at block boundaries too. CUDA uses production
    // defaults; also run with LOOKAHEAD=2 to exercise concurrent GPU BLAS.
    if (!std::getenv("MPLAPACK_LAPACK_LOOKAHEAD"))
        setenv("MPLAPACK_LAPACK_LOOKAHEAD", "2", 1);
#endif
#ifdef MPLAPACK_FACTOR_CUDA
    int devices=0;
    if (cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) return 77;
#endif
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    mpfr_set_default_prec(512);
#endif
    const char *env=std::getenv("MPLAPACK_LAPACK_FAST");
    fast=env && std::strcmp(env,"1")==0;
    INTEGER nb=iMlaenv(1,"Rgetrf"," ",600,600,-1,-1);
    std::set<INTEGER> sizes={1,2,3,nb-1,nb,nb+1,2*nb+3,257,599};
    for (INTEGER n:sizes) {
        if (n<1) continue;
        for (auto m: {n,n+7,std::max(INTEGER(1),n-7)}) {
            test_lu<REAL>(m,n,0); test_lu<COMPLEX>(m,n,0);
        }
        for (int kind=1;kind<=4;++kind) { test_lu<REAL>(n,n,kind); test_lu<COMPLEX>(n,n,kind); }
        for (const char *up:{"U","L"}) for (int bad=0;bad<2;++bad) {
            test_chol<REAL>(n,up,bad); test_chol<COMPLEX>(n,up,bad);
        }
    }
    // Tall, skinny panels reach the parallel row-block pivot reduction and
    // scaling, including ties spanning several blocks and a short tail.
    for (int kind=0; kind<=4; ++kind) {
        test_lu<REAL>(4103,3,kind);
        test_lu<COMPLEX>(4103,3,kind);
    }
    const bool saved_fast = fast;
    fast = false; // Public GETRF2 always retains the reference recursion.
    for (int kind=0; kind<=4; ++kind) {
        test_lu<REAL>(137,129,kind,true);
        test_lu<COMPLEX>(137,129,kind,true);
    }
    fast = saved_fast;
    test_signed_zeros();
#ifdef MPLAPACK_FACTOR_CUDA
    ++checks;
    cudaError_t status = cudaDeviceSynchronize();
    if (status != cudaSuccess) {
        ++failures;
        std::fprintf(stderr, "FAIL CUDA runtime (possible silent CPU fallback): %s\n", cudaGetErrorString(status));
    }
#endif
    std::printf("Rfactor: %d checks, %d failures, mode=%s, max backward error <= %s\n",
                checks,failures,fast?"fast":"default",worst<=REAL(30)?"30":"FAIL");
    return failures?1:0;
}
