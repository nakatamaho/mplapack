set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20"
set key font "Helvetica,16"
set title font "Helvetica,24"
set title "Rsyrk on Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz "
set xlabel "Dimension"
set ylabel "MFLOPS"
#set terminal postscript eps color enhanced
set terminal pdf

plot \
"log.Rsyrk._Float128"       using 1:3 title '\_Float128'             with lines linewidth 6, \
"log.Rsyrk._Float128_opt"   using 1:3 title '\_Float128 (OpenMP)'    with lines linewidth 6, \
"log.Rsyrk._Float64x"       using 1:3 title '\_Float64x'             with lines linewidth 6, \
"log.Rsyrk._Float64x_opt"   using 1:3 title '\_Float64x (OpenMP)'    with lines linewidth 6, \
"log.Rsyrk.dd"              using 1:3 title 'double double'          with lines linewidth 6, \
"log.Rsyrk.dd_opt"          using 1:3 title 'double double (OpenMP)' with lines linewidth 6, \
"log.Rsyrk.double"          using 1:3 title 'double'                 with lines linewidth 6

