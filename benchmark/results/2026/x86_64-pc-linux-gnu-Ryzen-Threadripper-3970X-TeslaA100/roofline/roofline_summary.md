| type | kernel | form | n | bound | % of roofline | % of GEMM plateau |
|---|---|---|---:|---|---:|---:|
| dd | cgemm | NC | 2048 | arith | 60.4 | 60.2 |
| dd | cgemm | NN | 2048 | arith | 60.7 | 60.5 |
| dd | cherk | LN | 1280 | arith | 1.4 | 1.4 |
| dd | cherk | LN_k64 | 2048 | arith | 1.8 | 1.8 |
| dd | ctrsm | LLNU | 1536 | arith | 3.1 | 3.1 |
| dd | ctrsm | RLCN_p64 | 2048 | arith | 5.6 | 5.6 |
| dd | gemm | NN | 2048 | arith | 99.8 | 99.4 |
| dd | gemm | NT | 2048 | arith | 99.8 | 99.4 |
| dd | gemm | TN | 2048 | arith | 100.1 | 99.8 |
| dd | gemm_sweep | NN | 2048 | arith | 99.8 |  |
| dd | gemm_update | NN | 2048 | arith | 95.5 | 95.2 |
| dd | getrf2 | P64 | 2048 | arith | 0.7 | 0.7 |
| dd | laswp | R64 | 2048 | memory | 5.5 |  |
| dd | lauum | L | 2048 | arith | 1.5 | 1.5 |
| dd | potrf2 | L | 256 | arith | 0.1 | 0.1 |
| dd | syrk | LN | 2048 | arith | 50.5 | 50.3 |
| dd | syrk | LN_k64 | 2048 | arith | 49.2 | 49.1 |
| dd | syrk | UT | 2048 | arith | 92.5 | 92.2 |
| dd | syrk | UT_k64 | 2048 | arith | 85.1 | 84.8 |
| dd | trsm | LLNU | 2048 | arith | 1.4 | 1.4 |
| dd | trsm | LLNU_p64 | 2048 | arith | 3.1 | 3.1 |
| dd | trsm | LUNN | 2048 | arith | 1.2 | 1.2 |
| dd | trsm | LUNN_p64 | 2048 | arith | 2.9 | 2.9 |
| dd | trsm | LUTN | 2048 | arith | 5.1 | 5.1 |
| dd | trsm | LUTN_p64 | 2048 | arith | 5.7 | 5.7 |
| dd | trsm | RLTN | 2048 | arith | 1.8 | 1.8 |
| dd | trsm | RLTN_p64 | 2048 | arith | 3.4 | 3.4 |
| dd | trtri | UN | 1280 | arith | 0.1 | 0.1 |
| td | cgemm | NC | 2048 | arith | 59.8 | 59.0 |
| td | cgemm | NN | 2048 | arith | 59.9 | 59.1 |
| td | cherk | LN | 1024 | arith | 2.2 | 2.2 |
| td | cherk | LN_k64 | 2048 | arith | 3.8 | 3.7 |
| td | ctrsm | LLNU | 1280 | arith | 4.8 | 4.8 |
| td | ctrsm | RLCN_p64 | 2048 | arith | 8.6 | 8.5 |
| td | gemm | NN | 2048 | arith | 101.0 | 99.6 |
| td | gemm | NT | 2048 | arith | 101.1 | 99.7 |
| td | gemm | TN | 2048 | arith | 101.0 | 99.7 |
| td | gemm_sweep | NN | 2048 | arith | 101.1 |  |
| td | gemm_update | NN | 2048 | arith | 99.6 | 98.3 |
| td | getrf2 | P64 | 2048 | arith | 1.8 | 1.7 |
| td | laswp | R64 | 2048 | memory | 7.2 |  |
| td | lauum | L | 1792 | arith | 1.5 | 1.4 |
| td | potrf2 | L | 256 | arith | 0.3 | 0.3 |
| td | syrk | LN | 2048 | arith | 49.7 | 49.0 |
| td | syrk | LN_k64 | 2048 | arith | 49.1 | 48.4 |
| td | syrk | UT | 2048 | arith | 98.0 | 96.7 |
| td | syrk | UT_k64 | 2048 | arith | 94.3 | 93.0 |
| td | trsm | LLNU | 2048 | arith | 4.5 | 4.4 |
| td | trsm | LLNU_p64 | 2048 | arith | 6.3 | 6.3 |
| td | trsm | LUNN | 2048 | arith | 4.5 | 4.4 |
| td | trsm | LUNN_p64 | 2048 | arith | 5.8 | 5.8 |
| td | trsm | LUTN | 2048 | arith | 5.6 | 5.5 |
| td | trsm | LUTN_p64 | 2048 | arith | 7.1 | 7.0 |
| td | trsm | RLTN | 1792 | arith | 3.1 | 3.0 |
| td | trsm | RLTN_p64 | 2048 | arith | 5.3 | 5.2 |
| td | trtri | UN | 1024 | arith | 0.2 | 0.2 |
| qd | cgemm | NC | 2048 | arith | 57.4 | 57.6 |
| qd | cgemm | NN | 2048 | arith | 57.4 | 57.6 |
| qd | cherk | LN | 928 | arith | 2.5 | 2.5 |
| qd | cherk | LN_k64 | 2048 | arith | 5.3 | 5.3 |
| qd | ctrsm | LLNU | 1024 | arith | 4.9 | 4.9 |
| qd | ctrsm | RLCN_p64 | 2048 | arith | 9.0 | 9.0 |
| qd | gemm | NN | 2048 | arith | 99.4 | 99.7 |
| qd | gemm | NT | 2048 | arith | 99.4 | 99.7 |
| qd | gemm | TN | 2048 | arith | 99.5 | 99.8 |
| qd | gemm_sweep | NN | 2048 | arith | 99.4 |  |
| qd | gemm_update | NN | 2048 | arith | 98.8 | 99.2 |
| qd | getrf2 | P64 | 2048 | arith | 3.3 | 3.3 |
| qd | laswp | R64 | 2048 | memory | 8.5 |  |
| qd | lauum | L | 1280 | arith | 1.1 | 1.1 |
| qd | potrf2 | L | 256 | arith | 0.3 | 0.3 |
| qd | syrk | LN | 2048 | arith | 52.9 | 53.1 |
| qd | syrk | LN_k64 | 2048 | arith | 52.6 | 52.7 |
| qd | syrk | UT | 2048 | arith | 96.2 | 96.5 |
| qd | syrk | UT_k64 | 2048 | arith | 94.7 | 95.0 |
| qd | trsm | LLNU | 2048 | arith | 7.3 | 7.4 |
| qd | trsm | LLNU_p64 | 2048 | arith | 7.7 | 7.7 |
| qd | trsm | LUNN | 1536 | arith | 5.5 | 5.5 |
| qd | trsm | LUNN_p64 | 2048 | arith | 7.0 | 7.0 |
| qd | trsm | LUTN | 1792 | arith | 6.7 | 6.7 |
| qd | trsm | LUTN_p64 | 2048 | arith | 9.3 | 9.3 |
| qd | trsm | RLTN | 1536 | arith | 4.2 | 4.2 |
| qd | trsm | RLTN_p64 | 2048 | arith | 6.4 | 6.4 |
| qd | trtri | UN | 928 | arith | 0.2 | 0.2 |
| mpfr256 | cgemm | NC | 1024 | arith | 18.8 | 41.3 |
| mpfr256 | cgemm | NN | 1024 | arith | 19.4 | 42.8 |
| mpfr256 | cherk | LN | 352 | arith | 0.2 | 0.5 |
| mpfr256 | cherk | LN_k64 | 864 | arith | 0.5 | 1.1 |
| mpfr256 | ctrsm | LLNU | 544 | arith | 0.8 | 1.7 |
| mpfr256 | ctrsm | RLCN_p64 | 1024 | arith | 0.2 | 0.5 |
| mpfr256 | gemm | NN | 2048 | arith | 45.5 | 100.3 |
| mpfr256 | gemm | NT | 2048 | arith | 45.7 | 100.7 |
| mpfr256 | gemm | TN | 2048 | arith | 43.6 | 96.1 |
| mpfr256 | gemm_sweep | NN | 2048 | arith | 45.2 |  |
| mpfr256 | gemm_update | NN | 2048 | arith | 41.0 | 90.3 |
| mpfr256 | getrf2 | P64 | 2048 | arith | 1.6 | 3.5 |
| mpfr256 | potrf2 | L | 256 | arith | 0.2 | 0.5 |
| mpfr256 | syrk | LN | 2048 | arith | 24.8 | 54.6 |
| mpfr256 | syrk | LN_k64 | 2048 | arith | 22.3 | 49.1 |
| mpfr256 | syrk | UT | 2048 | arith | 31.6 | 69.6 |
| mpfr256 | syrk | UT_k64 | 2048 | arith | 26.6 | 58.6 |
| mpfr256 | trsm | LLNU | 1536 | arith | 2.3 | 5.1 |
| mpfr256 | trsm | LLNU_p64 | 2048 | arith | 3.3 | 7.4 |
| mpfr256 | trsm | LUNN | 928 | arith | 1.4 | 3.2 |
| mpfr256 | trsm | LUNN_p64 | 2048 | arith | 2.9 | 6.3 |
| mpfr256 | trsm | LUTN | 992 | arith | 1.7 | 3.8 |
| mpfr256 | trsm | LUTN_p64 | 2048 | arith | 3.2 | 7.0 |
| mpfr256 | trsm | RLTN | 928 | arith | 1.5 | 3.3 |
| mpfr256 | trsm | RLTN_p64 | 2048 | arith | 3.1 | 6.8 |
| mpfr512 | cgemm | NC | 800 | arith | 8.3 | 24.0 |
| mpfr512 | cgemm | NN | 800 | arith | 8.5 | 24.5 |
| mpfr512 | cherk | LN | 256 | arith | 0.1 | 0.4 |
| mpfr512 | cherk | LN_k64 | 608 | arith | 0.3 | 0.9 |
| mpfr512 | ctrsm | LLNU | 352 | arith | 0.4 | 1.1 |
| mpfr512 | ctrsm | RLCN_p64 | 1024 | arith | 0.4 | 1.1 |
| mpfr512 | gemm | NN | 2048 | arith | 34.2 | 98.8 |
| mpfr512 | gemm | NT | 1536 | arith | 33.9 | 98.0 |
| mpfr512 | gemm | TN | 1536 | arith | 32.0 | 92.4 |
| mpfr512 | gemm_sweep | NN | 2048 | arith | 34.0 |  |
| mpfr512 | gemm_update | NN | 2048 | arith | 31.9 | 92.1 |
| mpfr512 | getrf2 | P64 | 2048 | arith | 1.8 | 5.3 |
| mpfr512 | potrf2 | L | 256 | arith | 0.2 | 0.5 |
| mpfr512 | syrk | LN | 2048 | arith | 13.4 | 38.6 |
| mpfr512 | syrk | LN_k64 | 2048 | arith | 14.0 | 40.5 |
| mpfr512 | syrk | UT | 2048 | arith | 20.5 | 59.2 |
| mpfr512 | syrk | UT_k64 | 2048 | arith | 16.6 | 48.0 |
| mpfr512 | trsm | LLNU | 1280 | arith | 2.0 | 5.7 |
| mpfr512 | trsm | LLNU_p64 | 2048 | arith | 3.3 | 9.5 |
| mpfr512 | trsm | LUNN | 736 | arith | 1.1 | 3.3 |
| mpfr512 | trsm | LUNN_p64 | 2048 | arith | 2.4 | 7.1 |
| mpfr512 | trsm | LUTN | 736 | arith | 1.3 | 3.7 |
| mpfr512 | trsm | LUTN_p64 | 2048 | arith | 2.6 | 7.5 |
| mpfr512 | trsm | RLTN | 736 | arith | 1.2 | 3.5 |
| mpfr512 | trsm | RLTN_p64 | 2048 | arith | 2.6 | 7.4 |
| mpfr1024 | cgemm | NC | 608 | arith | 8.5 | 59.6 |
| mpfr1024 | cgemm | NN | 544 | arith | 9.2 | 64.1 |
| mpfr1024 | cherk | LN | 160 | arith | 0.1 | 0.6 |
| mpfr1024 | cherk | LN_k64 | 416 | arith | 0.2 | 1.6 |
| mpfr1024 | ctrsm | LLNU | 224 | arith | 0.3 | 1.8 |
| mpfr1024 | ctrsm | RLCN_p64 | 1024 | arith | 0.6 | 3.9 |
| mpfr1024 | gemm | NN | 1024 | arith | 13.6 | 95.3 |
| mpfr1024 | gemm | NT | 1024 | arith | 13.6 | 95.4 |
| mpfr1024 | gemm | TN | 1024 | arith | 13.2 | 92.0 |
| mpfr1024 | gemm_sweep | NN | 1024 | arith | 13.6 |  |
| mpfr1024 | gemm_update | NN | 1024 | arith | 14.7 | 103.1 |
| mpfr1024 | getrf2 | P64 | 1024 | arith | 0.9 | 6.4 |
| mpfr1024 | potrf2 | L | 256 | arith | 0.1 | 0.8 |
| mpfr1024 | syrk | LN | 1024 | arith | 6.1 | 42.5 |
| mpfr1024 | syrk | LN_k64 | 1024 | arith | 7.4 | 51.5 |
| mpfr1024 | syrk | UT | 1024 | arith | 11.8 | 82.3 |
| mpfr1024 | syrk | UT_k64 | 1024 | arith | 9.9 | 69.5 |
| mpfr1024 | trsm | LLNU | 928 | arith | 1.5 | 10.8 |
| mpfr1024 | trsm | LLNU_p64 | 1024 | arith | 1.7 | 12.0 |
| mpfr1024 | trsm | LUNN | 544 | arith | 0.8 | 5.7 |
| mpfr1024 | trsm | LUNN_p64 | 1024 | arith | 1.1 | 7.7 |
| mpfr1024 | trsm | LUTN | 544 | arith | 0.9 | 6.3 |
| mpfr1024 | trsm | LUTN_p64 | 1024 | arith | 1.2 | 8.0 |
| mpfr1024 | trsm | RLTN | 544 | arith | 0.8 | 5.6 |
| mpfr1024 | trsm | RLTN_p64 | 1024 | arith | 1.1 | 7.5 |
| mpfr200rt | cgemm | NC | 64 | memory | - |  |
| mpfr200rt | cgemm | NN | 64 | memory | - |  |
| mpfr200rt | cherk | LN | 64 | memory | - |  |
| mpfr200rt | cherk | LN_k64 | 64 | memory | - |  |
| mpfr200rt | ctrsm | LLNU | 64 | memory | - |  |
| mpfr200rt | ctrsm | RLCN_p64 | 64 | memory | - |  |
| mpfr200rt | gemm | NN | 1024 | memory | - |  |
| mpfr200rt | gemm | NT | 864 | memory | - |  |
| mpfr200rt | gemm | TN | 864 | memory | - |  |
| mpfr200rt | gemm_update | NN | 1024 | memory | - |  |
| mpfr200rt | getrf2 | P64 | 1024 | memory | - |  |
| mpfr200rt | potrf2 | L | 256 | memory | - |  |
| mpfr200rt | syrk | LN | 1024 | memory | - |  |
| mpfr200rt | syrk | LN_k64 | 1024 | memory | - |  |
| mpfr200rt | syrk | UT | 928 | memory | - |  |
| mpfr200rt | syrk | UT_k64 | 1024 | memory | - |  |
| mpfr200rt | trsm | LLNU | 800 | memory | - |  |
| mpfr200rt | trsm | LLNU_p64 | 1024 | memory | - |  |
| mpfr200rt | trsm | LUNN | 416 | memory | - |  |
| mpfr200rt | trsm | LUNN_p64 | 1024 | memory | - |  |
| mpfr200rt | trsm | LUTN | 512 | memory | - |  |
| mpfr200rt | trsm | LUTN_p64 | 1024 | memory | - |  |
| mpfr200rt | trsm | RLTN | 416 | memory | - |  |
| mpfr200rt | trsm | RLTN_p64 | 1024 | memory | - |  |
