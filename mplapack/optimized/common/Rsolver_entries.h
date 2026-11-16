// Public wrappers shared by the four CUDA flavors. SOLVER_CALL is internal.
#ifndef MPLAPACK_SOLVER_COMPLEX_ONLY
bool SOLVER_ENTRY(Rgetrs)(const char *trans,INTEGER const n,INTEGER const nrhs,REAL *a,INTEGER const lda,INTEGER *ipiv,REAL *b,INTEGER const ldb,INTEGER &info) {
    return SOLVER_CALL("Rgetrs",false,false,trans,n,nrhs,a,lda,ipiv,b,ldb,info);
}
bool SOLVER_ENTRY(Rpotrs)(const char *uplo,INTEGER const n,INTEGER const nrhs,REAL *a,INTEGER const lda,REAL *b,INTEGER const ldb,INTEGER &info) {
    return SOLVER_CALL("Rpotrs",false,true,uplo,n,nrhs,a,lda,nullptr,b,ldb,info);
}
bool SOLVER_ENTRY(Rgesv)(INTEGER const n,INTEGER const nrhs,REAL *a,INTEGER const lda,INTEGER *ipiv,REAL *b,INTEGER const ldb,INTEGER &info) {
    return SOLVER_CALL("Rgesv",true,false,"N",n,nrhs,a,lda,ipiv,b,ldb,info);
}
bool SOLVER_ENTRY(Rposv)(const char *uplo,INTEGER const n,INTEGER const nrhs,REAL *a,INTEGER const lda,REAL *b,INTEGER const ldb,INTEGER &info) {
    return SOLVER_CALL("Rposv",true,true,uplo,n,nrhs,a,lda,nullptr,b,ldb,info);
}
#endif
#ifndef MPLAPACK_SOLVER_REAL_ONLY
bool SOLVER_ENTRY(Cgetrs)(const char *trans,INTEGER const n,INTEGER const nrhs,COMPLEX *a,INTEGER const lda,INTEGER *ipiv,COMPLEX *b,INTEGER const ldb,INTEGER &info) {
    return SOLVER_CALL("Cgetrs",false,false,trans,n,nrhs,a,lda,ipiv,b,ldb,info);
}
bool SOLVER_ENTRY(Cpotrs)(const char *uplo,INTEGER const n,INTEGER const nrhs,COMPLEX *a,INTEGER const lda,COMPLEX *b,INTEGER const ldb,INTEGER &info) {
    return SOLVER_CALL("Cpotrs",false,true,uplo,n,nrhs,a,lda,nullptr,b,ldb,info);
}
bool SOLVER_ENTRY(Cgesv)(INTEGER const n,INTEGER const nrhs,COMPLEX *a,INTEGER const lda,INTEGER *ipiv,COMPLEX *b,INTEGER const ldb,INTEGER &info) {
    return SOLVER_CALL("Cgesv",true,false,"N",n,nrhs,a,lda,ipiv,b,ldb,info);
}
bool SOLVER_ENTRY(Cposv)(const char *uplo,INTEGER const n,INTEGER const nrhs,COMPLEX *a,INTEGER const lda,COMPLEX *b,INTEGER const ldb,INTEGER &info) {
    return SOLVER_CALL("Cposv",true,true,uplo,n,nrhs,a,lda,nullptr,b,ldb,info);
}
#endif
#undef SOLVER_ENTRY
#undef SOLVER_CALL
