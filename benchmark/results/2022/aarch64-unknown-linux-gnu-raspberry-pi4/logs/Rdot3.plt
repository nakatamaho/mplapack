set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,16"
set title font "Helvetica,24"
set title "Rdot on "
set xlabel "Dimension"
set ylabel "MFLOPS"
#set terminal postscript eps color enhanced
set terminal pdf

plot "log.Rdot.double_opt"	    using 1:2 title 'double (OpenMP)'     with lines linewidth 6