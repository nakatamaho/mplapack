set title "Rlaruv test (0,1) uniform random"
set ylabel "count"
binwidth=0.01
bin(x,binwidth)=binwidth*floor(x/binwidth) + binwidth/2.0
set term png
set output 'Rlaruv.png'
plot 'Rlaruv.txt' using (bin($1,binwidth)):(1.0) smooth freq with boxes

