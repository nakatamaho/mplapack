    double *a_d = new double [ m * n];
    double *afac_d = new double [ m * n];
    double *work_d = new double [ lwork * 10] ;
    double *t2_d = new double[nb2 * n];
    double *q_d = new double[l * l];
    int lwork_d = (int)lwork;
    int info_d;
    for ( int p=1; p<=m; p++){
      for ( int q=1; q<= n; q++){
	a_d[ (p -1) + (q-1) * lda ] = 	a[ (p -1) + (q-1) * lda ];
	afac_d[ (p -1) + (q-1) * lda ] = a[ (p -1) + (q-1) * lda ];	
      }
    }
    strncpy(srnamt, "Rgetsqrhrt", srnamt_len);
    Rgetsqrhrt(m, n, mb1, nb1, nb2, af, m, t2, nb2, work, lwork, info);
    LAPACKE_dgetsqrhrt(LAPACK_COL_MAJOR, (int)m, (int)n, mb1, nb1, nb2, afac_d, m, t2_d, nb2);
    printf("a=");printmat(m,n,a,lda);printf("\n");
    printf("a_d=");printmat(m,n,a_d,lda);printf("\n");
    printf("af=");printmat(m,n,af,lda);printf("\n");
    printf("afac_d=");printmat(m,n,afac_d,lda);printf("\n");
    printf("t2=");printmat(nb2,n,t2,nb2);printf("\n");
    printf("t2_d=");printmat(nb2,n,t2_d,nb2);printf("\n");            
    printf("q=");printmat(m,m,q,ldq);printf("\n");

    LAPACKE_dgemqrt(LAPACK_COL_MAJOR, 'L', 'N', (int)m, (int)m, (int)k, (int)nb2_ub, afac_d, (int)m, t2_d,(int)nb2, q_d, (int)m);
    printf("q_d=");printmat(m,m,q_d,ldq);printf("\n");    

    delete []a_d;
    delete []afac_d;
    delete []work_d;
    delete []t2_d;
    delete []q_d;

