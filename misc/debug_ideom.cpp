        {
	  double *t_d = new double[ldt * n];
	  double *work_d = new double[max((INTEGER)1, lwork)];
	  double scale_d;
          double dtmp;
	  int ierr_d;
	  for (int pp=0; pp < ldt * n; pp++) t_d[pp] = t[pp];
	  for (int pp=0; pp < max(1,lwork); pp++) work_d[pp] = work[pp];
	  printf("# n1 %d n2 %d \n", (int)n1, (int)n2);
	  printf("a_d ="); printmat(n1, n1, t_d, ldt); printf("\n");
	  printf("b_d ="); printmat(n2, n2, &t_d[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], ldt); printf("\n");
	  printf("c_d ="); printmat(n1, n2, work_d, n1); printf("\n");	  	  
          LAPACKE_dtrsyl(LAPACK_COL_MAJOR, 'N', 'N', -1, (int)n1, (int)n2, t_d, (int)ldt, &t_d[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], (int)ldt, work_d, (int)n1, &scale_d);
	  printf("scale_d = "); printnum(scale_d);printf("\n");
	  printf("x_d ="); printmat(n1, n2, work_d, n1); printf("\n");
	  delete []t_d;
	  delete []work_d;
	}
        printf("a ="); printmat(n1, n1, t, ldt); printf("\n");
        printf("b ="); printmat(n2, n2, &t[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], ldt); printf("\n");
        printf("c ="); printmat(n1, n2, work, n1); printf("\n");
	Rtrsyl("N", "N", -1, n1, n2, t, ldt, &t[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], ldt, work, n1, scale, ierr);
        printf("scale = "); printnum_short(scale);printf("\n");	  
        printf("x ="); printmat(n1, n2, work, n1);printf("\n");

        {   
	  printf("dtrsen start\n");
	  double *t_d = new double[ldt * ldt];
	  double *q_d = new double[ldt * ldq];
	  double *wi_d = new double[n];
	  double *wr_d = new double[n];
	  double s_d;
	  double sep_d;
	  int select_d[ldt];
	  int m_d;
	  for (int pp = 1; pp <= n; pp++) {
	    for (int qq = 1; qq <= n; qq++) {
	      t_d[(pp - 1) + (qq - 1) * ldt] = t[(pp - 1) + (qq - 1) * ldt];
	      q_d[(pp - 1) + (qq - 1) * ldq] = q[(pp - 1) + (qq - 1) * ldq];
	    }
	  }
	  for (int pp = 1; pp <= n; pp++) {
	    select_d[pp-1]=select[pp-1];
	  }
	  LAPACKE_dtrsen(LAPACK_COL_MAJOR, 'B', 'V', select_d, (int)n, t_d, (int)ldt, q_d, (int)ldt, wr_d, wi_d, &m_d, &s_d, &sep_d);
	  printf("#n %d s_d, sep_d = ", count); printnum_short (s_d); printf(" "); printnum_short(sep_d); printf("\n");
	  for (int i = 1; i <= n; i = i + 1) {
            printf("dtrsen: wr_d, wi_d %d = ", (int)i);
            printnum(wr_d[i - 1]); printf(" "); printnum(wi_d[i - 1]);
            printf("\n");
	  }
	  delete []t_d; 
	  delete []q_d;
	  delete []wi_d;
	  delete []wr_d;
	  printf("dtrsen end\n");
	}	

        {
            double *t_d = new double[n * n];
            double *q_d = new double[n * n];
            double *wi_d = new double[n];
            double *wr_d = new double[n];
            for (int pp = 1; pp <= n; pp++) {
                for (int qq = 1; qq <= n; qq++) {
                    t_d[(pp - 1) + (qq - 1) * ldt] = t[(pp - 1) + (qq - 1) * ldt];
                    q_d[(pp - 1) + (qq - 1) * ldq] = q[(pp - 1) + (qq - 1) * ldq];
                }
            }
	    //            printf("tshur_d=");    printmat(n, n, t, ldt);     printf("\n");
            LAPACKE_dhseqr(LAPACK_COL_MAJOR, 'S', 'V', (int)n, 1, (int)n, t_d, (int)ldt, wr_d, wi_d, q_d, (int)ldq);
            for (int i = 1; i <= n; i = i + 1) {
                printf("wr_d, wi_d %d = ", (int)i);
                printnum(wr_d[i - 1]); printf(" "); printnum(wi_d[i - 1]);
                printf("\n");
            }
            delete[] wr_d;
            delete[] wi_d;
            delete[] t_d;
            delete[] q_d;
        }


{
   	    double *t_d = new double [ldt * n];
   	    double *q_d = new double [ldq * n];
	    double s_d, sep_d;
	    int m_d;
	    double *wrtmp_d = new double[n];
	    double *witmp_d = new double[n];
	    int *select_d = new int[n];	    
	    for ( int pp=1; pp<=n; pp++){
	      for ( int qq=1; qq<= n; qq++){
		t_d[ (pp -1) + (qq-1) * ldt ] = t[ (pp -1) + (qq-1) * ldt ];
		q_d[ (pp -1) + (qq-1) * ldq ] = q[ (pp -1) + (qq-1) * ldt ];		
	      }
	    }
  	    LAPACKE_dtrsen (LAPACK_COL_MAJOR, 'B', 'V', select_d, (int)n, t_d, (int)ldt, q_d, (int)ldt, wrtmp_d, witmp_d, &m_d, &s_d, &sep_d);
	    printf("sep_d = "); printnum(sep_d); printf("\n");
            delete [] select_d;
            delete [] witmp_d;
            delete [] wrtmp_d;
            delete [] q_d;
            delete [] t_d;
	}


/*
	{
          double *t_d = new double [n * n];
          double *work_d = new double [lwork];
          for ( int p=1; p<=n; p++){
            for ( int q=1; q<= n; q++){
              t_d[ (p -1) + (q-1) * ldt ] = t[ (p -1) + (q-1) * ldt ];
            }
          }
	  printf("torg=");  printmat(n, n, t, ldt);  printf("\n");
          LAPACKE_dgehrd(LAPACK_COL_MAJOR, (int)n, 1, n, t_d, (int)ldt, &work_d[1-1]);
          printf("tout_d=");printmat(n,n,t_d,ldt);printf("\n");
          delete [] t_d;
          delete [] work_d;
	}
	*/


    __complex__ double *afac_d;
    __complex__ double *tau_d;
    {
        afac_d = new __complex__ double[std::max(int(lda * n), 1)];
        tau_d = new __complex__ double[std::max(int(max(m, n)), 1)];
        for (int p = 1; p <= m; p++) {
            for (int q = 1; q <= n; q++) {
                __real__ afac_d[(p - 1) + (q - 1) * lda] = a[(p - 1) + (q - 1) * lda].real();
                __imag__ afac_d[(p - 1) + (q - 1) * lda] = a[(p - 1) + (q - 1) * lda].imag();
            }
        }
        LAPACKE_zgelqf(LAPACK_COL_MAJOR, m, n, afac_d, lda, tau_d);
        printf("af_d=");
        printmat(m, n, afac_d, lda);
        printf("\n");
    }
    delete[] afac_d;
    delete[] tau_d;


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
///////////////////////////////////////////////
                                        printf("ca=");   printnum(ca);              printf("\n");
                                        printf("a=");    printmat(na, na, a, lda);  printf("\n");
                                        printf("b=");    printmat(na, nw, b, ldb);  printf("\n");
                                        printf("w=");    printnum(wr);              printf("\n");
                                        printf("d1=");   printnum(d1);              printf("\n");
                                        printf("d2=");   printnum(d2);              printf("\n");
                                        printf("s=");    printnum(scale);           printf("\n");
                                        printf("X=");    printmat(na, 1, x, ldx);   printf("\n");
                                        printf("\n");
                                        printf("(ca * a - w * d) * X - (s * b)\n");
                                        printf("d=diag([d1,d2])\n");
                                        printf("info = %d\n", (int)info);
                                        bool trans = ltrans[itrans - 1];
                                        int na_d = (int)na;
                                        int nw_d = (int)nw;
                                        double smin_d = smin;
                                        double ca_d = ca;
                                        double a_d[4];
                                        a_d[0] = a[0];
                                        a_d[1] = a[1];
                                        a_d[2] = a[2];
                                        a_d[3] = a[3];
                                        int lda_d = 2;
                                        double d1_d = d1;
                                        double d2_d = d2;
                                        double b_d[4];
                                        b_d[0] = b[0];
                                        b_d[1] = b[1];
                                        b_d[2] = b[2];
                                        b_d[3] = b[3];
                                        double x_d[4];
                                        double wr_d = wr;
                                        double wi_d = wi;
                                        double scale_d;
                                        double xnorm_d;
                                        int info_d;
                                        dlaln2_(&trans, &na_d, &nw_d, &smin_d, &ca_d, a_d, &lda_d, &d1_d, &d2_d, b_d, &lda_d, &wr_d, &wi_d, x_d, &lda_d, &scale_d, &xnorm_d, &info_d);
                                        printf("ca_d=");   printnum(ca_d);              printf("\n");
                                        printf("a_d=");    printmat(na, na, a_d, lda);  printf("\n");
                                        printf("b_d=");    printmat(na, nw, b_d, ldb);  printf("\n");
                                        printf("w_d=");    printnum(wr_d);              printf("\n");
                                        printf("d1_d=");   printnum(d1_d);              printf("\n");
                                        printf("d2_d=");   printnum(d2_d);              printf("\n");
                                        printf("s_d=");    printnum(scale_d);           printf("\n");
                                        printf("X_d=");    printmat(na, 1, x_d, ldx);   printf("\n");
                                        printf("\n");
///////////////////////////////////////////////


	/*
	  double *t_d = new double [ldt * ldt];
	  double *le_d = new double [ldt * ldt];
	  double *re_d = new double [ldt * ldt];
	  int *select_d = new int [10];
          int m_d;
	  for (int p=1; p<=n; p++){
            for (int q=1; q<= n; q++) {
              t_d [ (p -1) + (q-1) * ldt ] = t[(p -1) + (q-1) * ldt];
            }
	  }
	  printf("t_d="); printmat(n, n, t, ldt); printf("\n");
	  LAPACKE_dtrevc(LAPACK_COL_MAJOR, 'B', 'A', select_d, (int)n, t_d, (int) ldt, le_d, ldt, re_d, (int)ldt, (int)n, &m_d);
	  printf("vlout_d="); printmat(n, n, le_d, ldt); printf("\n");
	  printf("vrout_d="); printmat(n, n, re_d, ldt); printf("\n");
	  delete [] select_d;
	  delete [] re_d;
	  delete [] le_d;
	  delete [] t_d;
	*/
