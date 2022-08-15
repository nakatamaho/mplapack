set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,16"
set title font "Helvetica,24"
set title "Rpotrf on Intel(R) Xeon(R) CPU E5-2623 v3 @ 3.00GHz "
set xlabel "Dimension"
set ylabel "MFLOPS"
#set terminal postscript eps color enhanced
set terminal pdf

plot \
"log.Rpotrf.double_opt"	    using 1:2 title 'double (OpenMP)'        with lines linewidth 6
