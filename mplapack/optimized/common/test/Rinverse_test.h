// Inverse regression extension of Rfactor_test; generated oracles unchanged.
void Rtrtri_reference_cpu(const char *,const char *,INTEGER,REAL *,INTEGER,INTEGER &);
void Ctrtri_reference_cpu(const char *,const char *,INTEGER,COMPLEX *,INTEGER,INTEGER &);
void Rlauum_reference_cpu(const char *,INTEGER,REAL *,INTEGER,INTEGER &);
void Clauum_reference_cpu(const char *,INTEGER,COMPLEX *,INTEGER,INTEGER &);
namespace {
#define Rtrtri Rtrtri_reference_cpu
#define Ctrtri Ctrtri_reference_cpu
#define Rlauum Rlauum_reference_cpu
#define Clauum Clauum_reference_cpu
#define Rpotri inverse_Rpotri_reference
#define Cpotri inverse_Cpotri_reference
#define Rgetri inverse_Rgetri_reference
#define Cgetri inverse_Cgetri_reference
#include "../../../reference/Rpotri.cpp"
#include "../../../reference/Cpotri.cpp"
#include "../../../reference/Rgetri.cpp"
#include "../../../reference/Cgetri.cpp"
#undef Rtrtri
#undef Ctrtri
#undef Rlauum
#undef Clauum
#undef Rpotri
#undef Cpotri
#undef Rgetri
#undef Cgetri
#define INVERSE_TEST_OVERLOAD(T,P) \
void inverse_tri(const char *u,const char *d,INTEGER n,T *a,INTEGER ld,INTEGER &i,bool ref) { if(ref) P##trtri_reference_cpu(u,d,n,a,ld,i); else P##trtri(u,d,n,a,ld,i); } \
void inverse_lauum(const char *u,INTEGER n,T *a,INTEGER ld,INTEGER &i,bool ref) { if(ref) P##lauum_reference_cpu(u,n,a,ld,i); else P##lauum(u,n,a,ld,i); } \
void inverse_potri(const char *u,INTEGER n,T *a,INTEGER ld,INTEGER &i,bool ref) { if(ref) inverse_##P##potri_reference(u,n,a,ld,i); else P##potri(u,n,a,ld,i); } \
void inverse_getri(INTEGER n,T *a,INTEGER ld,INTEGER *p,T *w,INTEGER lw,INTEGER &i,bool ref) { if(ref) inverse_##P##getri_reference(n,a,ld,p,w,lw,i); else P##getri(n,a,ld,p,w,lw,i); }
INVERSE_TEST_OVERLOAD(REAL,R)
INVERSE_TEST_OVERLOAD(COMPLEX,C)
#undef INVERSE_TEST_OVERLOAD

template<class T> std::vector<T> inverse_product(const std::vector<T> &a,const std::vector<T> &b,INTEGER n) {
    std::vector<T> c(size_t(n)*n,T(0));
    if(n) product(n,n,n,a.data(),b.data(),c.data());
    return c;
}
template<class T> std::vector<T> inverse_adjoint(const std::vector<T> &a,INTEGER n) {
    std::vector<T> b(a.size());
    for(INTEGER j=0;j<n;++j) for(INTEGER i=0;i<n;++i) b[i+j*n]=adjoint(a[j+i*n]);
    return b;
}
template<class T> void inverse_residual(const std::vector<T> &a,const std::vector<T> &b,INTEGER n) {
    auto p=inverse_product(a,b,n); REAL norm=0,err=0;
    for(INTEGER j=0;j<n;++j) {
        REAL col=0,ecol=0;
        for(INTEGER i=0;i<n;++i) { col+=magnitude(a[i+j*n]); ecol+=magnitude(p[i+j*n]-T(i==j?1:0)); }
        if(col>norm) norm=col; if(ecol>err) err=ecol;
    }
    REAL bnorm=0;
    for(INTEGER j=0;j<n;++j) { REAL col=0; for(INTEGER i=0;i<n;++i) col+=magnitude(b[i+j*n]); if(col>bnorm) bnorm=col; }
    REAL ratio=n?REAL(err/(REAL(n)*norm*bnorm*Rlamch("E"))):REAL(0);
    check(ratio<=REAL(30),"inverse residual",n,n,0);
}
template<class T> void inverse_case(INTEGER n,const char *up,bool unit,int kind) {
    const INTEGER ld=n+3; const bool upper=*up=='U';
    std::vector<T> a(size_t(ld)*n,T(19)),tri(size_t(n)*n,T(0));
    for(INTEGER j=0;j<n;++j) for(INTEGER i=0;i<n;++i) if(upper?i<=j:i>=j) {
        T v=i==j?T(2):T(random_scalar((T*)nullptr)/T(std::max(INTEGER(1),n)));
        if(kind==2 && i!=j) v=T(-0.0);
        a[i+j*ld]=v; tri[i+j*n]=i==j && unit?T(1):v;
    }
    if(kind==1 && n) a[n/2+(n/2)*ld]=T(0);
    auto r=a,g=a; INTEGER ir=-1,ig=-2;
    routine=std::is_same<T,REAL>::value?"Rtrtri":"Ctrtri"; triangle=*up;
    inverse_tri(up,unit?"U":"N",n,r.data(),ld,ir,true);
    inverse_tri(up,unit?"U":"N",n,g.data(),ld,ig,false);
    check(ir==ig,"info",n,n,kind);
    if(!fast) identical(r,g,n,n,kind);
    if(!ig && n) {
        std::vector<T> inv(size_t(n)*n,T(0));
        for(INTEGER j=0;j<n;++j) for(INTEGER i=0;i<n;++i) if(upper?i<=j:i>=j) inv[i+j*n]=i==j && unit?T(1):g[i+j*ld];
        inverse_residual(tri,inv,n);
    }
    if(unit || kind==1) return;
    routine=std::is_same<T,REAL>::value?"Rlauum":"Clauum";
    r=a;g=a;inverse_lauum(up,n,r.data(),ld,ir,true);inverse_lauum(up,n,g.data(),ld,ig,false);
    check(ir==ig,"info",n,n,kind); if(!fast) identical(r,g,n,n,kind);
    if(fast && n) {
        auto t=inverse_adjoint(tri,n),expected=upper?inverse_product(tri,t,n):inverse_product(t,tri,n);
        std::vector<T> actual(size_t(n)*n,T(0));
        for(INTEGER j=0;j<n;++j) for(INTEGER i=0;i<n;++i) actual[i+j*n]=(upper?i<=j:i>=j)?g[i+j*ld]:adjoint(g[j+i*ld]);
        residual(expected,actual,n,n,n,kind);
    }
    routine=std::is_same<T,REAL>::value?"Rpotri":"Cpotri";
    r=a;g=a;inverse_potri(up,n,r.data(),ld,ir,true);inverse_potri(up,n,g.data(),ld,ig,false);
    check(ir==ig,"info",n,n,kind);if(!fast) identical(r,g,n,n,kind);
    if(!ig && n) {
        auto t=inverse_adjoint(tri,n),original=upper?inverse_product(t,tri,n):inverse_product(tri,t,n);
        std::vector<T> inv(size_t(n)*n,T(0));
        for(INTEGER j=0;j<n;++j) for(INTEGER i=0;i<n;++i) inv[i+j*n]=(upper?i<=j:i>=j)?g[i+j*ld]:adjoint(g[j+i*ld]);
        inverse_residual(original,inv,n);
    }
}
template<class T> void inverse_general(INTEGER n,bool singular) {
    routine=std::is_same<T,REAL>::value?"Rgetri":"Cgetri";triangle='-';
    INTEGER ld=n+3,nb=iMlaenv(1,routine," ",n,-1,-1,-1),lw=std::max(INTEGER(1),n*nb);
    std::vector<T> a(size_t(ld)*n,T(19)),l(size_t(n)*n,T(0)),u(l.size(),T(0));
    std::vector<INTEGER> p(n,n);
    for(INTEGER j=0;j<n;++j) for(INTEGER i=0;i<n;++i) {
        T v=i==j?T(2):T(random_scalar((T*)nullptr)/T(std::max(INTEGER(1),n)));
        a[i+j*ld]=v;if(i>j) l[i+j*n]=v; else u[i+j*n]=v;
        if(i==j) l[i+j*n]=T(1);
    }
    if(singular && n) a[n/2+(n/2)*ld]=T(0);
    auto r=a,g=a;std::vector<T> wr(lw,T(17)),wg=wr;INTEGER ir=-1,ig=-2;
    inverse_getri(n,r.data(),ld,p.data(),wr.data(),-1,ir,true);
    inverse_getri(n,g.data(),ld,p.data(),wg.data(),-1,ig,false);
    check(ir==ig && same(wr[0],wg[0]),"workspace query",n,n,0);identical(a,g,n,n,0);
    inverse_getri(n,r.data(),ld,p.data(),wr.data(),lw,ir,true);
    inverse_getri(n,g.data(),ld,p.data(),wg.data(),lw,ig,false);
    check(ir==ig && same(wr[0],wg[0]),"info/workspace",n,n,singular);
    if(!fast) { identical(r,g,n,n,singular); identical(wr,wg,n,n,singular); }
    if(!ig && n) {
        auto original=inverse_product(l,u,n);std::vector<T> inv(size_t(n)*n);
        for(INTEGER i=n;i-->0;) for(INTEGER j=0;j<n;++j) std::swap(original[i+j*n],original[p[i]-1+j*n]);
        for(INTEGER j=0;j<n;++j) for(INTEGER i=0;i<n;++i) inv[i+j*n]=g[i+j*ld];
        inverse_residual(original,inv,n);
    }
}
void test_inverses(bool quick) {
    std::vector<INTEGER> sizes=quick?std::vector<INTEGER>{0,3,65}:std::vector<INTEGER>{0,1,2,3,31,32,33,63,64,65,131,257,599};
    for(INTEGER n:sizes) {
        for(const char *up:{"U","L"}) for(bool unit:{false,true}) for(int kind:{0,1,2}) {
            inverse_case<REAL>(n,up,unit,kind);inverse_case<COMPLEX>(n,up,unit,kind);
        }
        for(bool singular:{false,true}) { inverse_general<REAL>(n,singular);inverse_general<COMPLEX>(n,singular); }
    }
}
}
