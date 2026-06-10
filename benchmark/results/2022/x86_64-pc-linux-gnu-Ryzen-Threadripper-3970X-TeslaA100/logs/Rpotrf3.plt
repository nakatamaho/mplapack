set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,16"
set title font "Helvetica,20"
set title "Rpotrf on AMD Ryzen Threadripper 3970X 32-Core Processor "
set xlabel "Dimension"
set ylabel "MFLOPS"
set key above
set terminal pdf

plot \
"log.Rpotrf.double_opt"	    using 1:2 title 'double (OpenMP)'        with lines linewidth 1, \
"log.Rpotrf.double"         using 1:2 title 'double'                 with lines linewidth 1, \
"log.dpotrf.ref"	    using 1:2 title 'double (Ref.BLAS)'      with lines linewidth 1, \
"log.dpotrf.openblas"       using 1:2 title 'double (OpenBLAS)'      with lines linewidth 1
