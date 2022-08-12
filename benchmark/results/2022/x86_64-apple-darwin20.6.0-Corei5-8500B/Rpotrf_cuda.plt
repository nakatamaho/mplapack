set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,20"
set title font "Helvetica,24"
set title "Rpotrf performance on GPU"
set xlabel "Dimension"
set ylabel "MFLOPS"
set terminal postscript eps color enhanced
set yrange [0:40000]

plot "log.Rpotrf.dd_cuda_total"  using 1:4 title 'Total ' with lines linewidth 1
