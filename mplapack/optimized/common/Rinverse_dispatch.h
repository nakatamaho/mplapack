// Included by basename-shadowing BF inverse dispatchers. CPU fallback is the
// exact shared optimized CPU implementation, with CPU-only BLAS and leaves.
#include <mpblas.h>
#include <mplapack.h>
#if defined(MPLAPACK_BUILD_WITH_DD)
#include <mpblas_dd_cuda.h>
#define INVERSE_PUBLIC_ENTRY_IMPL(R) R##_dd_cuda
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include <mpblas_td_cuda.h>
#define INVERSE_PUBLIC_ENTRY_IMPL(R) R##_td_cuda
#else
#include <mpblas_qd_cuda.h>
#define INVERSE_PUBLIC_ENTRY_IMPL(R) R##_qd_cuda
#endif
#define INVERSE_PUBLIC_ENTRY(R) INVERSE_PUBLIC_ENTRY_IMPL(R)
#if MPLAPACK_INVERSE_COMPLEX
using InverseHost=COMPLEX;
#else
using InverseHost=REAL;
#endif
#define MPLAPACK_TRTRI_CPU_DISPATCH
#define MPLAPACK_LAUUM_CPU_DISPATCH
#if MPLAPACK_INVERSE_KIND == 0
#if MPLAPACK_INVERSE_COMPLEX
#include "Ctrtri.cpp"
#define INVERSE_PUBLIC Ctrtri
#else
#include "Rtrtri.cpp"
#define INVERSE_PUBLIC Rtrtri
#endif
void INVERSE_PUBLIC(const char *u,const char *d,INTEGER const n,InverseHost *a,INTEGER const lda,INTEGER &info) {
    if(!INVERSE_PUBLIC_ENTRY(INVERSE_PUBLIC)(u,d,n,a,lda,info))trtri_cpu(u,d,n,a,lda,info);
}
#elif MPLAPACK_INVERSE_KIND == 1
#if MPLAPACK_INVERSE_COMPLEX
#include "Clauum.cpp"
#define INVERSE_PUBLIC Clauum
#else
#include "Rlauum.cpp"
#define INVERSE_PUBLIC Rlauum
#endif
void INVERSE_PUBLIC(const char *u,INTEGER const n,InverseHost *a,INTEGER const lda,INTEGER &info) {
    if(!INVERSE_PUBLIC_ENTRY(INVERSE_PUBLIC)(u,n,a,lda,info))lauum_cpu(u,n,a,lda,info);
}
#else
// GETRI/POTRI retain the generated outer algorithm, but bind its inverse
// calls directly to the same optimized CPU engine, not another GPU attempt.
#include "Rinverse_cpu_blas.h"
#include "Rinverse_cpu_common.h"
namespace {
#if MPLAPACK_INVERSE_COMPLEX
#include "Ctrtri.cpp"
#if MPLAPACK_INVERSE_KIND == 2
#include "Clauum.cpp"
#define Ctrtri trtri_cpu
#define Clauum lauum_cpu
#define Cpotri potri_cpu
#include "../../reference/Cpotri.cpp"
#undef Ctrtri
#undef Clauum
#undef Cpotri
#define INVERSE_PUBLIC Cpotri
#else
#define Ctrtri trtri_cpu
#define Cgetri getri_cpu
#include "../../reference/Cgetri.cpp"
#undef Ctrtri
#undef Cgetri
#define INVERSE_PUBLIC Cgetri
#endif
#else
#include "Rtrtri.cpp"
#if MPLAPACK_INVERSE_KIND == 2
#include "Rlauum.cpp"
#define Rtrtri trtri_cpu
#define Rlauum lauum_cpu
#define Rpotri potri_cpu
#include "../../reference/Rpotri.cpp"
#undef Rtrtri
#undef Rlauum
#undef Rpotri
#define INVERSE_PUBLIC Rpotri
#else
#define Rtrtri trtri_cpu
#define Rgetri getri_cpu
#include "../../reference/Rgetri.cpp"
#undef Rtrtri
#undef Rgetri
#define INVERSE_PUBLIC Rgetri
#endif
#endif
}
#if MPLAPACK_INVERSE_KIND == 2
void INVERSE_PUBLIC(const char *u,INTEGER const n,InverseHost *a,INTEGER const lda,INTEGER &info) {
    if(!INVERSE_PUBLIC_ENTRY(INVERSE_PUBLIC)(u,n,a,lda,info))potri_cpu(u,n,a,lda,info);
}
#else
void INVERSE_PUBLIC(INTEGER const n,InverseHost *a,INTEGER const lda,INTEGER *p,InverseHost *work,INTEGER const lwork,INTEGER &info) {
    if(!INVERSE_PUBLIC_ENTRY(INVERSE_PUBLIC)(n,a,lda,p,work,lwork,info))getri_cpu(n,a,lda,p,work,lwork,info);
}
#endif
#endif
#undef INVERSE_PUBLIC
#undef INVERSE_PUBLIC_ENTRY
#undef INVERSE_PUBLIC_ENTRY_IMPL
#undef MPLAPACK_TRTRI_CPU_DISPATCH
#undef MPLAPACK_LAUUM_CPU_DISPATCH
