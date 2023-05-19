set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,16"
set title font "Helvetica,20"
set title "Rgemv on AMD Ryzen Threadripper 3970X 32-Core Processor "
set xlabel "Dimension"
set ylabel "MFLOPS"
set key above
set terminal pdf

plot \
"log.Rgemv.double_opt"	    using 1:3 title 'double (OpenMP)'     with lines linewidth 1, \
"log.Rgemv.double"          using 1:3 title 'double'              with lines linewidth 1, \
"log.dgemv.ref"             using 1:3 title 'double (Ref.BLAS)'   with lines linewidth 1, \
"log.dgemv.openblas"        using 1:3 title 'double (OpenBLAS)'   with lines linewidth 1
