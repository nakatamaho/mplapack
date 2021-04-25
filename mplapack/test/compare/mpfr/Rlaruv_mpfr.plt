set encoding utf8
set xlabel "i-th random number"
set ylabel "(0,1)
set term png
set output 'Rlaruv_mpfr.png'
plot  "Rlaruv_mpfr.txt" using 1:2
