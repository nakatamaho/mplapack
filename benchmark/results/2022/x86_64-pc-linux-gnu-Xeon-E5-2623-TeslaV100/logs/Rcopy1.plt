set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,16"
set title font "Helvetica,24"
set title "Rcopy on Intel(R) Xeon(R) CPU E5-2623 v3 @ 3.00GHz "
set xlabel "Dimension"
set ylabel "MFLOPS"
set terminal pdf

plot \
"log.Rcopy._Float128"       using 1:2 title '\_Float128'             ,\
"log.Rcopy._Float128_opt"   using 1:2 title '\_Float128 (OpenMP)'    ,\
"log.Rcopy._Float64x"       using 1:2 title '\_Float64x'             ,\
"log.Rcopy._Float64x_opt"   using 1:2 title '\_Float64x (OpenMP)'    ,\
"log.Rcopy.dd"              using 1:2 title 'double double'          ,\
"log.Rcopy.dd_opt"          using 1:2 title 'double double (OpenMP)' ,\
"log.Rcopy.double"          using 1:2 title 'double'                 

