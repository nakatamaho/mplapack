set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,16"
set title font "Helvetica,20"
set title "Rgetrf on %%MODELNAME%%"
set xlabel "Dimension"
set ylabel "MFLOPS"
#set terminal postscript eps color enhanced
set terminal pdf

plot \
"log.Rgetrf.double_opt"	    using 1:3 title 'double (OpenMP)'        with lines linewidth 1
