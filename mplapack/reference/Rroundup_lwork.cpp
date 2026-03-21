REAL Rroundup_lwork(INTEGER const lwork) {
    REAL rlwork = castREAL(lwork);
    if (castINTEGER(rlwork) < lwork) {
        rlwork = rlwork * (1.0 + Rlamch("E"));
    }
    return rlwork;
}
