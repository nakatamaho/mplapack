// Included once by each backend bridge after defining blas1_entry<Complex>.
// Public declarations come exclusively from the generated accelerator header.
#define BLAS1_REAL_ENTRY(NAME,OP) \
bool MPLAPACK_CUDA_ENTRY(NAME)(INTEGER const n,REAL *x,INTEGER const incx,REAL &result) { \
    REAL r;COMPLEX c;long index=0; \
    if(!blas1_entry<false>(mplapack_blas1::OP,n,x,incx,(REAL*)nullptr,1,r,c,index)) return false; \
    result=r;return true; }
#define BLAS1_COMPLEX_REAL_ENTRY(NAME,OP) \
bool MPLAPACK_CUDA_ENTRY(NAME)(INTEGER const n,COMPLEX *x,INTEGER const incx,REAL &result) { \
    REAL r;COMPLEX c;long index=0; \
    if(!blas1_entry<true>(mplapack_blas1::OP,n,x,incx,(COMPLEX*)nullptr,1,r,c,index)) return false; \
    result=r;return true; }
BLAS1_REAL_ENTRY(Rasum,asum)
BLAS1_REAL_ENTRY(Rnrm2,norm)
BLAS1_COMPLEX_REAL_ENTRY(RCasum,asum)
BLAS1_COMPLEX_REAL_ENTRY(RCnrm2,norm)
#undef BLAS1_REAL_ENTRY
#undef BLAS1_COMPLEX_REAL_ENTRY
bool MPLAPACK_CUDA_ENTRY(Rdot)(INTEGER const n,REAL *x,INTEGER const incx,REAL *y,INTEGER const incy,REAL &result) {
    REAL r;COMPLEX c;long index=0;
    if(!blas1_entry<false>(mplapack_blas1::dot,n,x,incx,y,incy,r,c,index)) return false;
    result=r;return true;
}
#define BLAS1_COMPLEX_DOT(NAME,OP) \
bool MPLAPACK_CUDA_ENTRY(NAME)(INTEGER const n,COMPLEX *x,INTEGER const incx,COMPLEX *y,INTEGER const incy,COMPLEX &result) { \
    REAL r;COMPLEX c;long index=0; \
    if(!blas1_entry<true>(mplapack_blas1::OP,n,x,incx,y,incy,r,c,index)) return false; \
    result=c;return true; }
BLAS1_COMPLEX_DOT(Cdotc,dotc)
BLAS1_COMPLEX_DOT(Cdotu,dotu)
#undef BLAS1_COMPLEX_DOT
#define BLAS1_IAMAX(NAME,TYPE,ISCOMPLEX) \
bool MPLAPACK_CUDA_ENTRY(NAME)(INTEGER const n,TYPE *x,INTEGER const incx,INTEGER &result) { \
    REAL r;COMPLEX c;long index=0; \
    if(!blas1_entry<ISCOMPLEX>(mplapack_blas1::iamax,n,x,incx,(TYPE*)nullptr,1,r,c,index)) return false; \
    result=index+1;return true; }
BLAS1_IAMAX(iRamax,REAL,false)
BLAS1_IAMAX(iCamax,COMPLEX,true)
#undef BLAS1_IAMAX
