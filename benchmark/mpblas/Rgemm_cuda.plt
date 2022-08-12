set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,20"
set title font "Helvetica,24"
set title "Rgemm GPU"
set xlabel "Dimension"
set ylabel "MFLOPS"
#set terminal postscript eps color enhanced
set terminal pdf

plot "log.Rgemm.dd_cuda_kernel" using 1:4 title 'Kernel only' with lines linewidth 1,\
     "log.Rgemm.dd_cuda_total"  using 1:4 title 'Total ' with lines linewidth 1
