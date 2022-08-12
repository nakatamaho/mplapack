set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,20"
set title font "Helvetica,24"
set title "Rsyrk performance on NVIDIA Tesla V100 PCIe 32GB] (rev a1)"
set xlabel "Dimension"
set ylabel "MFLOPS"
#set terminal postscript eps color enhanced
set terminal pdf

plot "log.Rsyrk.dd_cuda_kernel" using 1:3 title 'Kernel only' with lines linewidth 1,\
     "log.Rsyrk.dd_cuda_total"  using 1:3 title 'Total ' with lines linewidth 1
