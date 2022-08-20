set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,16"
set title font "Helvetica,20"
set title "Rpotrf on Tesla V100-PCIE-32GB"
set xlabel "Dimension"
set ylabel "MFLOPS"
#set terminal postscript eps color enhanced
set terminal pdf

plot "log.Rpotrf.dd_cuda_total"  using 1:2 title 'Total' with lines linewidth 1
