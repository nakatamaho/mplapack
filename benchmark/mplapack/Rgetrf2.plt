set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,16"
set title font "Helvetica,24"
set title "Rgetrf on %%MODELNAME%%"
set xlabel "Dimension"
set ylabel "MFLOPS"
#set terminal postscript eps color enhanced
set terminal pdf

plot \
"log.Rgetrf.mpfr"            using 1:2 title 'MPFR 512bit'            with lines linewidth 6, \
"log.Rgetrf.mpfr_opt"        using 1:2 title 'MPFR 512bit(OpenMP)'    with lines linewidth 6, \
"log.Rgetrf.gmp"             using 1:2 title 'GMP  512bit'            with lines linewidth 6, \
"log.Rgetrf.gmp_opt"         using 1:2 title 'GMP  512bit(OpenMP)'    with lines linewidth 6, \
"log.Rgetrf.qd"              using 1:2 title 'quad-double'            with lines linewidth 6, \
"log.Rgetrf.qd_opt"	    using 1:2 title 'quad-double(OpenMP)'    with lines linewidth 6
