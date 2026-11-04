// Included after the private tri_entry implementation and generated API header.
#define MPLAPACK_TRI_ENTRIES(P,H,C) \
bool MPLAPACK_CUDA_ENTRY(P##trsv)(const char *u,const char *t,const char *d,INTEGER const n,H *a,INTEGER const lda,H *x,INTEGER const inc) { return tri_entry<C>(0,u,t,d,n,0,a,lda,x,inc); } \
bool MPLAPACK_CUDA_ENTRY(P##tbsv)(const char *u,const char *t,const char *d,INTEGER const n,INTEGER const k,H *a,INTEGER const lda,H *x,INTEGER const inc) { return tri_entry<C>(1,u,t,d,n,k,a,lda,x,inc); } \
bool MPLAPACK_CUDA_ENTRY(P##tpsv)(const char *u,const char *t,const char *d,INTEGER const n,H *a,H *x,INTEGER const inc) { return tri_entry<C>(2,u,t,d,n,0,a,0,x,inc); }
MPLAPACK_TRI_ENTRIES(R,REAL,false)
MPLAPACK_TRI_ENTRIES(C,COMPLEX,true)
#undef MPLAPACK_TRI_ENTRIES
