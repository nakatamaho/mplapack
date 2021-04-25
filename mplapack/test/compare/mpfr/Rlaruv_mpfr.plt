set encoding utf8
set xlabel "i-th random number"
set ylabel "(0,1)
set title "Rlaruv test (0,1) uniform random"
set term png
set output 'Rlaruv_mpfr.png'
plot  "Rlaruv_mpfr.txt"
